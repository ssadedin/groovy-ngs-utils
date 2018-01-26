package gngs

import groovy.transform.CompileStatic
import htsjdk.samtools.SAMFileReader
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
class OrderedPairReader {
    
    SAMFileReader randomLookupReader
    
    List<SAMRecordPair> writeSpool = new LinkedList()
    Set<String> preprocessedReads = new HashSet()
    Map<String,SAMRecordPair> readIndex
    
    int pairs = 0
        
    boolean includeUnmapped 
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
    
    OrderedPairReader(Map options=[:], SAM sam) {
        this.randomLookupReader = sam.newReader()
            
        this.includeUnmapped = ((boolean)options?.includeUnmapped)
        this.sample = sam.getSamples()[0] ?: "";
            
        // The spool holds a buffer that tries to pair up reads before outputting them
        // The right size for it depends on the coverage depth and the separation between
        // paired reads. It should be able to hold enough reads that most of the time
        // a read's pair is encountered before the buffer overflows. When it overflows
        // the loop below starts doing random queries (slow) to resolve the mates for reads
        // so that the pair can be output together. 
        this.maxSpoolSize = options.spoolSize ? (int)options.spoolSize : 16000i
        this.bam = sam
        this.verbose = false
        this.readIndex = new HashMap(this.maxSpoolSize*2)
        
        this.progress = new ProgressCounter(
            extra:{
                "$sample ${currentRegion?:currentChr}, Spooled=${writeSpool.size()} (${writeSpool[-1].minPos()}-${writeSpool[0].minPos()}), Preproc=${preprocessedReads.size()} Unpaired=${unpairedReads}, ForcedQueries=$forcedQueries, Pairs=${pairs}"
            }, withRate: true, timeInterval:10000)
    }
    
    void eachPair(Closure c) {
        
       SAMFileReader r = bam.newReader() 
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
            randomLookupReader.close()
            iter.close()
        } 
    }
    
    @CompileStatic
    void eachPair(Regions regions, Closure c) {
        
        if(regions == null) {
            eachPair(c)
            return
        }
        
        SAMFileReader reader = this.bam.newReader()
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
            randomLookupReader.close()
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
                        
            SAMRecordPair pair = readIndex[readName]
            if(pair == null) {
                if(verbose)
                    println "XX: Buffer read $readName $record.referenceName:$record.alignmentStart (spool position ${writeSpool.size()})"
                pair = new SAMRecordPair(r1:record)
                writeSpool << pair
                readIndex[readName] = pair
                continue                    
            }
                
            if(verbose)
                println "XX: Pair found for $readName ($record.referenceName:$record.alignmentStart)"
                    
            pair.r2 = record
                        
            if(writeSpool.size()>maxSpoolSize) {
                forceQueryOverflowRead(record.alignmentStart)
            }
       
            processPartneredReads(c)
        }
    
        List<SAMRecordPair> spoolResidue = new LinkedList()
        for(SAMRecordPair pair in writeSpool) {
            if(pair.hasBothReads() && pair.r1.referenceName == currentChr) {
                c(pair.r1, pair.r2)
                readIndex.remove(pair.r1.readName)
                ++pairs                        
            }
            else {
                spoolResidue << pair
            }
        }
        writeSpool = spoolResidue
    }
    
    boolean isProcessableRead(SAMRecord record) {
//        println "Proper pair = " + record.getProperPairFlag()
//        println "Read paired = " + record.readPairedFlag
//                            
        
        String readName = record.readName
        
        // Unmapped reads are not proper pairs, so for those we do not apply the
        // check about proper pairs. For 
        if(includeUnmapped && (record.readUnmappedFlag || record.mateUnmappedFlag)) {
            if(!record.getReadPairedFlag())
                return false
        }
        else
        if(!record.getReadPairedFlag() || record.getReadUnmappedFlag() || !record.getProperPairFlag())
            return false
                            
        if((record.referenceIndex != record.mateReferenceIndex) || record.isSecondaryOrSupplementary()) {
            if(verbose)
                println "Chimeric / non-primary alignment $readName :" + record.referenceName + ":" + record.alignmentStart
            return false
        }
                                
        if(readName in preprocessedReads) {
//                        println "XX: Already processed $r1.readName"
            preprocessedReads.remove(readName)
            return false
        }
        
        return true
    }
    
    void processPartneredReads(Closure c) {
        while(writeSpool[0]?.hasBothReads()) {
            SAMRecordPair pair = writeSpool.remove(0)
            String writeReadName = pair.r1.readName
            if(verbose)
                println "XX: write pair $writeReadName ($pair.r1.referenceName:$pair.r1.alignmentStart, $pair.r2.referenceName:$pair.r2.alignmentStart)"
            c(pair.r1, pair.r2)
            ++pairs
            readIndex.remove(writeReadName)
        }
    }
    
    void forceQueryOverflowRead(int position) {
        
        SAMRecordPair pair = writeSpool[0]
        
        if(pair.r2 != null)
            return
        
        SAMRecord mate = null
        
        // If our position is far past the known position of the read's mate, assume it is NOT in the BAM file 
        // and leave it null
        if(position < pair.r1.mateAlignmentStart + 1000) {
            mate = bam.queryMate(randomLookupReader,pair.r1) 
            ++forcedQueries
        }
        else {
            ++unpairedReads
        }
            
        if(mate != null) {
            if(verbose)
                println "XX: Pair found for ${writeSpool[0].r1.readName} (${writeSpool[0].r1.alignmentStart}, ${mate.alignmentStart}) by random lookup "
            writeSpool[0].r2 = mate
            preprocessedReads.add(mate.readName)
        }
        else {
            SAMRecord mateless = writeSpool.remove(0).r1
            readIndex.remove(mateless.readName)
            if(verbose)
                println "XX: Read $mateless.readName has no mates :-("
        }
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
        readIndex.clear()
    }
}
