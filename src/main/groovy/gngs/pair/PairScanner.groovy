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
package gngs.pair

import gngs.*

import static gngs.Utils.human

import java.util.concurrent.ConcurrentHashMap

import groovy.transform.CompileStatic
import groovy.util.logging.Log
import groovyx.gpars.actor.Actor
import groovyx.gpars.actor.DefaultActor
import groovyx.gpars.group.DefaultPGroup
import groovyx.gpars.group.PGroup
import htsjdk.samtools.SamReader
import htsjdk.samtools.BAMIndex
import htsjdk.samtools.BAMIndexMetaData
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SAMRecordIterator
import htsjdk.samtools.SAMSequenceRecord

/**
 * Scans a BAM file and feeds reads encountered to a pool of {@link PairLocator}
 * instances for matching to pairs. The {@link PairLocator} is chosen based on
 * read name so that any given locator is guaranteed to see both a read and it's
 * mate, if the mate exists.
 * <p>
 * This class also supports a "sharding" parameter ({@link #shardId},
 * {@link #shardSize}, which allows only every nth read to be emitted with an
 * offset. The result is that you can run multiple instances of this class in
 * parallel and guarantee that each instance will emit a distinct set of reads.
 * 
 * @author Simon Sadedin
 */
@Log
class PairScanner {
    
    static PairScanner running = null
    
    SAM bam
    
    SAMRecord lastRead
    
    ProgressCounter progress = new ProgressCounter(withRate:true, withTime:true, timeInterval: 60000, extra: { 
        lastRead?.referenceName + ':' + lastRead?.alignmentStart + ", loc: " +  
        [locators*.received.sum(),locators*.paired.sum(),locators*.buffer*.size().sum()].collect {human(it)}.join(',') +
        " chm: ${human(locators*.chimeric.sum())}" + 
        " shf: ${human(shufflers[0].currentBufferSize)}" +
        " fmt: ${human(formatters*.formatted.sum())}, wrote: " + human(pairWriter.written) + 
        ", mem: ${Utils.human(CompactReadPair.memoryUsage)}/${Utils.human(Runtime.runtime.freeMemory())}/${Utils.human(Runtime.runtime.totalMemory())}"
    })
    
    PairWriter pairWriter
    
    PairWriter pairWriter2
    
    List<PairLocator> locatorIndex = []
    
    List<PairLocator> locators = []
    
    List<PairFilter> filters = []
    
    List<PairFormatter> formatters
    
    Regions regions
    
    List<Actor> actors = []
    
    int chimeric
    
    int shardId = -1
    
    int shardSize = 0
    
    int numLocators 
    
    boolean throttleWarning = false
    
    String filterExpr
    
    String debugRead = null
    
    Set<Integer> chromosomesWithReads = Collections.newSetFromMap(new ConcurrentHashMap(500))
    
    /**
     * How many formatted blocks ready to write will be buffered for writing.
     * <p>
     * Each block is roughly the size of the {@link #FORMATTER_BUFFER_SIZE} parameter
     */
    final static int DEFAULT_WRITER_QUEUE_SIZE = 1000
    
    int numFormatters = 2
    
    /**
     * The size of blocks to format before writing.
     * <p>
     * Each block is accumulated until it reaches thi size. This is essentially similar to 
     * output file buffering, but occurs before the file system layer.
     */
    final static int FORMATTER_BUFFER_SIZE = 1000 * 1000
    
    /**
     * The buffer within which to shuffle reads so as to randomise their output order
     */
    int shuffleBufferSize = 5000
    
    /**
     * The tag from which to extract base quality scores (use actual base qualities if null)
     */
    String baseQualityTag = null
    
    List<Shuffler> shufflers
    
    PGroup formatterGroup
    PGroup shufflerPGroup
    
    PairScanner(Writer writer1, Writer writer2, int numLocators, Regions regions = null, String filterExpr = null, int writerQueueSize = DEFAULT_WRITER_QUEUE_SIZE) {
        this.pairWriter = new PairWriter(writer1, writerQueueSize)
        this.pairWriter.parallelGroup = new DefaultPGroup(1)
        
        this.pairWriter2 = new PairWriter(writer2, writerQueueSize)
        
        this.regions = regions
        this.numLocators = numLocators
        this.filterExpr = filterExpr
        progress.log = log
    }
  
    PairScanner(Writer writer, int numLocators, Regions regions = null, String filterExpr = null, int writerQueueSize = DEFAULT_WRITER_QUEUE_SIZE) {
        this.pairWriter = new PairWriter(writer, writerQueueSize)
        this.regions = regions
        this.numLocators = numLocators
        this.filterExpr = filterExpr
        progress.log = log
    }
  
    void initLocators(SAM bam) {
        
        if(this.debugRead)
            formatters*.debugRead = this.debugRead
            
        int locatorsCreated = 0
        this.locators = []
        
        this.formatters = (1..numFormatters).collect { new PairFormatter(FORMATTER_BUFFER_SIZE, pairWriter, pairWriter2)  }
        
        this.formatterGroup = new DefaultPGroup(numFormatters+1)
        this.formatters*.parallelGroup = this.formatterGroup
        
        this.shufflers = this.formatters.collect { new Shuffler(it, shuffleBufferSize) }
        this.shufflerPGroup = new DefaultPGroup(this.shufflers.size()+1)
        this.shufflers*.parallelGroup = this.shufflerPGroup
        
        Set<Integer> sequencesWithReads = getContigsWithReads(bam)
        
        log.info "The following contigs have at least one read: " + sequencesWithReads.join(', ')
        
        // Note: since we may be running in sharded mode, certain locator positions
        // may never get used. We counter this by checking at each position
        // if the locator will be used and we keep creating locators until we 
        // have the requisite number at each of the shard positions in our array
        for(int i=0; locatorsCreated < numLocators; ++i) {
            if(shardId<0 || ((i%shardSize) == shardId)) {
                ++locatorsCreated
                
                createLocator(bam, sequencesWithReads, shufflers[i%shufflers.size()])
            }
            else {
                this.locatorIndex.add(null)
            }
        }
        
        // Fill up any trailing null positions needed
        int requiredIndexSize = shardSize*numLocators
        while(this.locatorIndex.size()<requiredIndexSize) {
            this.locatorIndex.add(null)
        }
        
        log.info "Created ${locatorsCreated} read pair locators"
        
        this.actors = locators + filters + shufflers + formatters + [
            pairWriter,
        ]
        
        if(pairWriter2 != null)
            this.actors << pairWriter2
            
        this.actors*.start()
    }
    
    @CompileStatic
    void createLocator(SAM bam, Set<Integer> sequencesWithReads, Shuffler shuffler) {
        
        PairLocator pl 
        if(filterExpr != null) {
            PairFilter filter = new PairFilter(shuffler, filterExpr) 
            this.filters << filter
            pl = new PairLocator(filter, sequencesWithReads)
            pl.compact = false // pass full SAMRecordPair through
        }
        else {
            pl = new PairLocator(shuffler, sequencesWithReads)
        }
        pl.baseQualityTag = this.baseQualityTag
                
        if(this.regions)
            pl.regions = new Regions((Iterable)this.regions)        
            
        if(debugRead != null)
            pl.debugRead = debugRead
            
        this.locators << pl
        this.locatorIndex << pl
    }
    
    /**
     * Interrogate the BAM index to determine which contigs have reads.
     * <p>
     * This is done to help better cope with BAM files where selected contigs have been 
     * included, leaving large numbers of mateless reads. By knowing up front that there
     * are no reads for a given contig, and therefore a mate positioned in that contig will
     * never be encountered, we can avoid storing those reads in memory.
     * 
     * @param bam   BAM file to check
     * @return  set of the indices of the reference sequences (contigs / chromosomes) that have 
     *          at least one read
     */
    Set<Integer> getContigsWithReads(SAM bam) {
        Set<Integer> sequencesWithReads = Collections.newSetFromMap(new ConcurrentHashMap())
        bam.withReader { SamReader r -> 
            List<Integer> seqIndices = r.fileHeader.sequenceDictionary.sequences*.sequenceIndex
            BAMIndex index = r.index
            for(Integer ind : seqIndices) {
                BAMIndexMetaData meta = index.getMetaData(ind)
                if(meta.getAlignedRecordCount() +  meta.getUnalignedRecordCount() + meta.getNoCoordinateRecordCount()>0) {
                    sequencesWithReads.add(ind)
                }
            }
        }
        
        return sequencesWithReads
    }
    
    @CompileStatic
    void scan(SAM bam) {
        log.info "Beginning scan of $bam.samFile"
        if(debugRead != null)
            log.info "Debugging read $debugRead"
            
        running = this
        this.initLocators(bam)
        try {
            this.scanBAM(bam)
        }
        finally {
            log.info "Stopping parallel threads ..."
            
            locators.eachWithIndex { a, i -> stopActor("Locator $i", a) }
            
            filters.eachWithIndex { a, i -> stopActor("Filter $i", a) }
            
            shufflers.eachWithIndex { a, i -> stopActor("Shuffler $i", a) }
            
            formatters.eachWithIndex { f, i ->
                stopActor "Formatter $i", f
            }
            stopActor "Writer", pairWriter
            if(pairWriter2)
                stopActor "Writer2", pairWriter2
                
            progress.end()
            running = null
        }
    }
    
    /**
     * Scan all the reads in the given BAM file
     */
    @CompileStatic
    private void scanBAM(SAM bam) {
        
        final SamReader reader = bam.newReader(fast:true)
        try {
            final SAMRecordIterator i = reader.iterator()
            try{
                scanBAMIterator(i)
            }
            finally {
                i.close()
            }
        }
        finally {
            reader.close()
        }        
    }

    /**
     * Process all reads yielded by the given iterator
     * <p>
     * Note: Caller must close iterator
     */
    @CompileStatic
    private void scanBAMIterator(final SAMRecordIterator i) {
        
        // These are in the hope that we get some compiler
        // optimisations - have not tested though
        final int locatorSize = locatorIndex.size()
        final int shardId = this.shardId
        final int shardSize = this.shardSize
        int readCount = 0
        while(i.hasNext()) {
            final SAMRecord read = i.next()
            lastRead = read
            progress.count()
            final int hash = Math.abs(read.readName.hashCode())
            final int locatorOffset = hash % locatorSize
            PairLocator locator = locatorIndex[locatorOffset]
            if(locator != null) {
                locator.assign(read)
            }
            else {
                if(!debugRead.is(null) && debugRead == read.readName)
                    log.info "Read $debugRead not assigned to shard (hash=$hash, lcoffset=$locatorOffset/$locatorSize)"
            }
            ++readCount
        }
        
        log.info "Scanned ${readCount} reads"
        
    }
    
    void stopActor(String name, Actor actor) {
        log.info "Stopping $name"
        if(actor instanceof RegulatingActor) {
            actor.sendStop()
        }
        else {
            actor << "stop"
        }
        actor.join()
    }
}
