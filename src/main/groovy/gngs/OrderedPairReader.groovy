/*
 *  Groovy NGS Utils - Some simple utilites for processing Next Generation Sequencing data.
 *
 *  Copyright (C) 2018 Simon Sadedin, ssadedin<at>gmail.com
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
package gngs

import java.io.IOException
import java.util.concurrent.atomic.AtomicInteger
import java.util.logging.Logger
import groovy.transform.CompileStatic
import groovy.util.logging.Log
import htsjdk.samtools.SamReader
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SAMRecordIterator

/**
 * A utility class for efficiently iterating all the read pairs in a BAM file.
 * <p>
 * NOTE: unpaired and chimeric reads are omitted
 * <p>
 * Performance depends on buffering a window of reads so that the mates for each read are
 * encountered within the buffer and thus random lookup of mates is not needed.  The right size 
 * for it depends on the coverage depth and the separation between
 * paired reads. It should be able to hold enough reads that most of the time
 * a read's pair is encountered before the buffer overflows. When it overflows
 * the loop below starts doing random queries (slow) to resolve the mates for reads
 * so that the pair can be output together. See the <code>spoolSize</code> option
 * which can be passed to the constructor.
 * 
 * @author Simon Sadedin
 */
@CompileStatic
class OrderedPairReader implements Closeable {
    
    SamReader randomLookupReader
    
    Logger log = Logger.getLogger("OrderedPairReader")
    
    /**
     * List of reads in order that need to be written out
     */
    List<SAMRecordPair> writeSpool = new LinkedList()
    
    Set<String> preprocessedReads = Collections.synchronizedSet(new HashSet())
    
    /**
     * To contain memory consumption, preprocessed reads are cleared every chromosome.
     * We can't, however, do that for chimeric reads since they span chromosomes
     */
    Set<String> chimericPreprocessedReads = Collections.synchronizedSet(new HashSet())
    
    MissingMateIndex missingMates
    
    int pairs = 0
        
    boolean includeUnmapped = false
    
    boolean includeChimeric = false
    
    String sample 
        
    // The spool holds a buffer that tries to pair up reads before outputting them
    // The right size for it depends on the coverage depth and the separation between
    // paired reads. It should be able to hold enough reads that most of the time
    // a read's pair is encountered before the buffer overflows. When it overflows
    // the loop below starts doing random queries (slow) to resolve the mates for reads
    // so that the pair can be output together. 
    final int maxSpoolSize 
        
    int forcedQueries = 0
    
    int unpairedReads = 0
    
    String currentChr = null 
    
    SAM bam
    
    boolean verbose
    
    ProgressCounter progress
    
    Region currentRegion = null
    
    int lookupConcurrency = 4
    int lookupReaderIndex = 0
    
    String debugRead = "H233NCCXX150123:6:2216:5912:67112"
    
    AtomicInteger asyncRequests = new AtomicInteger(0)
    
    AtomicInteger asyncResolves = new AtomicInteger(0)
    
    List<ReadMateLookupActor> lookupActors = []
    
    List<AsyncMateResolveWorker> resolvers = []
    
    ReadMateLookupActor forcedQueryLookupActor = null
    
    OrderedPairReader(Map options=[:], SAM sam) {
        this.randomLookupReader = sam.newReader()
        this.missingMates = new MissingMateIndex(this.maxSpoolSize*2)
        
        if('lookupConcurrency' in options) {
            this.lookupConcurrency = (int)options.lookupConcurrency
        }
        
        this.lookupActors = (1..lookupConcurrency).collect { 
            new ReadMateLookupActor(sam.newReader(), missingMates, { SAMRecordPair pair ->
                if(!pair.hasBothReads()) {
                    pair.flags.unpaired = true
                }
                this.addPreprocessedRead(pair)
                asyncResolves.getAndIncrement()
            })
        }

        this.resolvers = lookupActors.collect { 
            new AsyncMateResolveWorker(index:missingMates, lookupActor: it)
        }
        
        this.resolvers.eachWithIndex { r, index ->
            Thread t = new Thread(r, "Resolver-${index}")
            t.daemon = true
            t.start()
        }
        
        log.info "Started ${resolvers.size()} asynchronous resolvers"
        
        this.forcedQueryLookupActor = new ReadMateLookupActor(sam.newReader(), missingMates, { SAMRecordPair pair ->
            this.addPreprocessedRead(pair)
        })
        
//        this.lookupActors*.start()
        
        
        this.includeUnmapped = ((boolean)options?.includeUnmapped)
        this.sample = sam.getSamples()[0] ?: "";
            
        // The spool holds a buffer that tries to pair up reads before outputting them
        // The right size for it depends on the coverage depth and the separation between
        // paired reads. It should be able to hold enough reads that most of the time
        // a read's pair is encountered before the buffer overflows. When it overflows
        // the loop below starts doing random queries (slow) to resolve the mates for reads
        // so that the pair can be output together. 
        this.maxSpoolSize = options.spoolSize ? (int)options.spoolSize : 32000i
        this.bam = sam
        this.verbose = false
        this.progress = new ProgressCounter(
            extra:{
                
                String pos = writeSpool[0]?.r1Pos
                
//                "$sample ${currentRegion?:currentChr}, Spooled=${writeSpool.size()} (${writeSpool[-1].minPos()}-${writeSpool[0].minPos()}), Preproc=${preprocessedReads.size()} Lookaheads=${lookAheads},${lookAheadResolutions} Unpaired=${unpairedReads}, ForcedQueries=$forcedQueries, Pairs=${pairs}"
                "$sample $pos\t${currentRegion?:currentChr}, Spooled=${writeSpool.size()}, ${spoolStats()} Preproc=${preprocessedReads.size()} " + 
                "resolveQueue=${missingMates.resolveQueueSize}, " + 
                "Lookaheads=${forcedQueryLookupActor.lookAheads},${forcedQueryLookupActor.lookAheadResolutions} Unpaired=${unpairedReads}, ChimericBuffer=${missingMates.chimericReadIndex.size()} Forced=$forcedQueries, Async=${asyncResolves}/${asyncRequests} Pairs=${pairs}"
            }, withRate: true, timeInterval:10000)
    }
    
    void eachPair(Closure c) {
        
       SamReader r = bam.newReader() 
       try {
           eachPair(r.iterator(),c)
       }
       finally {
           r.close()
       }
    }
    
    void eachPair(SAMRecordIterator iter, Closure c) {
        
        try {
            this.eachPairImpl(iter,c)
        }
        finally {
            progress.end()
            iter.close()
            this.close()
        } 
    }
    
    @CompileStatic
    void eachPair(Regions regions, Closure c) {
        
        if(regions == null) {
            eachPair(c)
            return
        }
        
        SamReader reader = this.bam.newReader()
        try {
            for(Region region in regions) {
                currentRegion = region
                SAMRecordIterator<SAMRecord> iter = reader.query(region.chr, region.from, region.to, false)
                try {
                    eachPairImpl(iter,c)
                    
                    flushSpooledReads(c)
                }
                finally {
                    iter.close()
                }
            }
        }
        finally {
            progress.end()
            reader.close()
        }
    }
    
    void eachPairImpl(SAMRecordIterator iter, Closure c) {
        
        while(iter.hasNext()) {
            SAMRecord record = (SAMRecord)iter.next();
            String readName = record.readName
            if(isNewChromosome(record)) {
                flushSpooledReads(c)
            }
            progress.count()
                        
            if(!isProcessableRead(record))
                continue
                        
            currentChr = record.referenceName
            SAMRecordPair pair = missingMates.getPair(readName)
            if(pair == null) {
                if(verbose)
                    println "XX: Buffer read $readName $record.referenceName:$record.alignmentStart (spool position ${writeSpool.size()})"
                pair = new SAMRecordPair(r1:record)
                writeSpool << pair
                int matePos = record.mateAlignmentStart
                
                missingMates.add(pair)
                
//                // Chimeric reads are not going to get resolved on this chromosome, so
//                // feed them straight into asynchronous resolution
//                if(pair.chimeric && includeChimeric) {
//                    resolveAsynchronously(pair)
//                }
//                
            }
            else {
                    
                if(verbose)
                    println "XX: Pair found for $readName ($record.referenceName:$record.alignmentStart)"
                        
                pair.r2 = record
                missingMates.resolve(record.alignmentStart,pair)
            }
                
            if(writeSpool.size()>maxSpoolSize) {
                scheduleAsyncResolves()
                    
                forceQueryOverflowRead(record.alignmentStart)
            }
           
            processPartneredReads(c)
        }
    
        List<SAMRecordPair> spoolResidue = new LinkedList()
        for(SAMRecordPair pair in writeSpool) {
            if(pair.hasBothReads() && pair.r1.referenceName == currentChr) {
                c(pair.r1, pair.r2)
                ++pairs                        
            }
            else {
                spoolResidue << pair
            }
        }
        writeSpool = spoolResidue
    }
    
    @CompileStatic
    Map spoolStats() {
        writeSpool.take(100).countBy { SAMRecordPair pair ->
            if(pair.hasBothReads())
                "resolved"
            else
            if(pair.flags.unpaired) {
                "unpaired"
            }
            else
            if(pair.flags.requested) {
                "requested"
            }
            else {
                "pending"
            }
                
        }
    }
    
    void scheduleAsyncResolves() {
        // Let's async resolve some reads too
        int requested = 0
        int i = 1
        int prerequested = 0
        int unpaired = 0
        int presolved = 0
        for(; i<writeSpool.size() && i<100 && requested<20;) {
                    
            SAMRecordPair asyncPair = this.writeSpool[i++]
                    
            if(asyncPair.hasBothReads()) {
                ++presolved
                continue
            }
                    
            if(asyncPair.flags.unpaired) {
                ++unpaired
                continue
            }
                    
            if(asyncPair.flags.requested) {
                ++presolved
                continue
            }
                        
            if(this.resolveAsynchronously(asyncPair))
                ++requested
                        
        }
        if(verbose && requested>0)
            log.info "Requested ${requested} async resolves (prereq=$prerequested, unpaired=$unpaired, presolved=$presolved, index = $i/${writeSpool.size()}), queue size = ${missingMates.resolveQueue.size()}"
    }
    
    boolean resolveAsynchronously(SAMRecordPair pair) {
//        int actorIndex = lookupReaderIndex++ % lookupActors.size()
//        this.lookupActors[actorIndex] << new ResolveRequest(pair:pair)
        
        synchronized(this.missingMates) {
            this.missingMates.resolveQueue.add(pair)
            this.missingMates.notify()
        }
        pair.flags.requested = true
        this.asyncRequests.incrementAndGet()
        return true
    }
    
    boolean isProcessableRead(SAMRecord record) {
        
        if((debugRead !=null) && record.readName == debugRead) {
          log.info "XXX: Proper pair = " + record.getProperPairFlag()
          log.info "XXX: Read paired = " + record.readPairedFlag
          log.info "Include chimeric: " + this.includeChimeric
          log.info "Different mate index? " + (record.referenceIndex != record.mateReferenceIndex)
        }
        
        String readName = record.readName
        
        // Unmapped reads are not proper pairs, so for those we do not apply the
        // check about proper pairs. 
        if(includeUnmapped && (record.readUnmappedFlag || record.mateUnmappedFlag)) {
            if(!record.getReadPairedFlag())
                return false
        }
        else
        if(!record.getReadPairedFlag() || record.readUnmappedFlag || (!record.properPairFlag && !includeChimeric))
            return false
                            
        if((!includeChimeric && (record.referenceIndex != record.mateReferenceIndex)) || record.isSecondaryOrSupplementary()) {
            if(verbose)
                println "Chimeric / non-primary alignment $readName :" + record.referenceName + ":" + record.alignmentStart
                
            if((debugRead !=null) && record.readName == debugRead) {
                log.info "Chimeric / non-primary alignment $readName :" + record.referenceName + ":" + record.alignmentStart
            }
                
            return false
        }
        
                                
        if(readName in preprocessedReads) {
//                        println "XX: Already processed $r1.readName"
            preprocessedReads.remove(readName)
            return false
        }
        else
        if(readName in chimericPreprocessedReads) {
            chimericPreprocessedReads.remove(readName)
            return false
        }
        
        if((debugRead !=null) && record.readName == debugRead) {
           log.info "Including read $debugRead" 
        }
        
        return true
    }
    
    @CompileStatic
    void processPartneredReads(Closure c) {
        while(writeSpool[0]?.hasBothReads()) {
            SAMRecordPair pair = writeSpool.remove(0)
            assert pair.r2  != null
            String writeReadName = pair.r1.readName
            if(verbose)
                println "XX: write pair $writeReadName ($pair.r1.referenceName:$pair.r1.alignmentStart, $pair.r2.referenceName:$pair.r2.alignmentStart)"
            assert pair.r2  != null
            c(pair.r1, pair.r2)
            assert pair.r2  != null
            ++pairs
            assert pair.r2  != null
            assert missingMates != null
            missingMates.resolve(pair.r2.alignmentStart, pair)
        }
    }
    
    void forceQueryOverflowRead(int position) {
        
        SAMRecordPair pair = writeSpool[0]
        
        if(pair.r2 != null)
            return
            
        SAMRecord mate = null
        
        // If our position is far past the known position of the read's mate, assume it is NOT in the BAM file 
        // and leave it null
        int mateAlignmentStart = pair.r1.mateAlignmentStart
        if(!pair.flags.unpaired && position < mateAlignmentStart + 5000) {
            
            int regionScanSize = 500
            int matePosition = mateAlignmentStart
            
            RegionScan scanRegion = missingMates.allocateScannableRegion(pair, regionScanSize, 50)
            if(scanRegion != null) {
                mate = forcedQueryLookupActor.queryMateByRegionScan(pair.r1, scanRegion.region.region)
                scanRegion.markUnpaired()
            }
            else {
                mate = queryMateByDirectLookup(pair.r1)
            }
            
            ++forcedQueries
        }
        else {
            ++unpairedReads
        }
            
        missingMates.resolve(mateAlignmentStart, pair)
        
        if(mate != null) {
            if(verbose)
                println "XX: Pair found for ${writeSpool[0].r1.readName} (${writeSpool[0].r1.alignmentStart}, ${mate.alignmentStart}) by random lookup "
            writeSpool[0].r2 = mate
            
            addPreprocessedRead(writeSpool[0])
        }
        else {
            SAMRecord mateless = writeSpool.remove(0).r1
            if(verbose)
                println "XX: Read $mateless.readName has no mates :-("
        }
    }
    
    SAMRecord queryMateByDirectLookup(SAMRecord record) {
        bam.queryMate(randomLookupReader,record) 
    }
    
    boolean isNewChromosome(SAMRecord record) {
        currentChr != null && currentChr != record.referenceName        
    }
    
    void flushSpooledReads(Closure c) {
        for(SAMRecordPair readPair in writeSpool) {
            
            SAMRecord r1 = readPair.r1
            if(readPair.hasBothReads()) {
                c(readPair.r1, readPair.r2)
                progress.count()
            }
            else {
                ++unpairedReads
            }
        }
                            
        if(verbose)
            println "XX: ${unpairedReads} reads ignored because mate never encountered"
            
        writeSpool.clear()
        missingMates.clear() 
    }
    
    void addPreprocessedRead(SAMRecordPair pair) {
        if(pair.chimeric)
            chimericPreprocessedReads.add(pair.readName)
        else
            preprocessedReads.add(pair.readName)
    }

    @Override
    public void close() throws IOException {
        randomLookupReader.close()
        for(ReadMateLookupActor lookupActor in lookupActors) {
//            lookupActor << "stop"
//            lookupActor.join()
            lookupActor.close()
        }
    }
}
