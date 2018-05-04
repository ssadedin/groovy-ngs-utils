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


import static Utils.human

import groovy.transform.CompileStatic
import groovy.util.logging.Log
import groovyx.gpars.actor.Actor
import groovyx.gpars.actor.DefaultActor
import htsjdk.samtools.SAMFileReader
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SAMRecordIterator

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
    
    ProgressCounter progress = new ProgressCounter(withRate:true, withTime:true, extra: { 
        lastRead?.referenceName + ':' + lastRead?.alignmentStart + ", loc: " +  
        [locators*.received.sum(),locators*.paired.sum(),locators*.buffer*.size().sum()].collect {human(it)}.join(',') +
        " chimeric: ${human(locators*.chimeric.sum())}" + 
        " formatted: ${human(formatter.formatted)}, written: " + human(pairWriter.written)
    })
    
    PairWriter pairWriter
    
    List<PairLocator> locatorIndex = []
    
    List<PairLocator> locators = []
    
    List<PairFilter> filters = []
    
    PairLocator chimericLocator
    
    PairFormatter formatter 
    
    Regions regions
    
    List<Actor> actors = []
    
    int chimeric
    
    int shardId = -1
    
    int shardSize = 0
    
    int numLocators 
    
    boolean throttleWarning = false
    
    String filterExpr
    
    String debugRead = null
    
    /**
     * If more than this number of reads are unwritten, assume that we are
     * limited by downstream ability to consume our reads and start backing off
     */
    int maxWriteBufferSize = 3000000
    
    PairScanner(Writer writer, int numLocators, Regions regions = null, String filterExpr = null) {
        this.pairWriter = new PairWriter(writer)
        this.formatter = new PairFormatter(1000 * 1000, pairWriter)
        this.chimericLocator = new PairLocator(formatter)
        this.regions = regions
        this.numLocators = numLocators
        this.filterExpr = filterExpr
        progress.log = log
    }
    
    void initLocators() {
        
        if(this.debugRead)
            formatter.debugRead = this.debugRead
        
        int locatorsCreated = 0
        this.locators = []
        
        // Note: since we may be running in sharded mode, certain locator positions
        // may never get used. We counter this by checking at each position
        // if the locator will be used and we keep creating locators until we 
        // have the requisite number at each of the shard positions in our array
        for(int i=0; locatorsCreated < numLocators; ++i) {
            if(shardId<0 || ((i%shardSize) == shardId)) {
                ++locatorsCreated
                
                createLocator()
            }
            else {
                this.locatorIndex.add(null)
            }
        }
        
        this.actors = locators + filters + [
            chimericLocator,
            formatter,
            pairWriter,
        ]
        this.actors*.start()
    }
    
    @CompileStatic
    void createLocator() {
        PairLocator pl 
        if(filterExpr != null) {
            PairFilter filter = new PairFilter(formatter, filterExpr)
            this.filters << filter
            pl = new PairLocator(filter)
            pl.compact = false // pass full SAMRecordPair through
        }
        else {
            pl = new PairLocator(formatter)
        }
                
        if(this.regions)
            pl.regions = new Regions((Iterable)this.regions)        
            
        if(debugRead != null)
            pl.debugRead = debugRead
            
        this.locators << pl
        this.locatorIndex << pl
    }
    
    @CompileStatic
    void scan(SAM bam) {
        log.info "Beginning scan of $bam.samFile"
        if(debugRead != null)
            log.info "Debugging read $debugRead"
            
        running = this
        this.initLocators()
        try {
            this.scanBAM(bam)
        }
        finally {
            log.info "Stopping parallel threads ..."
            
            locators.eachWithIndex { a, i -> stopActor("Locator $i", a) }
            
            filters.eachWithIndex { a, i -> stopActor("Filter $i", a) }
            
            stopActor "Chimeric Locator", chimericLocator
            stopActor "Formatter", formatter
            stopActor "Writer", pairWriter
            progress.end()
            running = null
        }
    }
    
    private void scanBAM(SAM bam) {
        
        // These are in the hope that we get some compiler
        // optimisations - have not tested though
        final int locatorSize = locatorIndex.size()
        final int shardId = this.shardId
        final int shardSize = this.shardSize
        final int maxBufferedReads = this.maxWriteBufferSize
        
        final SAMFileReader reader = bam.newReader()
        reader.enableCrcChecking(false)
        reader.enableIndexCaching(true)
        reader.enableIndexMemoryMapping(true)
        
        try {
            final SAMRecordIterator i = reader.iterator()
            try {
                while(i.hasNext()) {
                    SAMRecord read = i.next()
                    lastRead = read
                    progress.count()
                    int hash = read.readName.hashCode()
                    int locatorOffset = hash % locatorSize
                    PairLocator locator = locatorIndex[locatorOffset]
                    if(locator != null) {
                        locator << read 
                        if(pairWriter.pending.get() > maxBufferedReads) {
                            if(!throttleWarning) {
                                log.info "Throttling output due to slow downstream consumption of reads"
                                throttleWarning = true
                            }
                            Thread.sleep(50)
                        }
                    }
                    else {
                        if(debugRead == read.readName)
                            log.info "Read $debugRead not assigned to shard (hash=$hash, lcoffset=$locatorOffset/$locatorSize)"
                    }
                        
                }
            }
            finally {
                i.close()
            }
        }
        finally {
            reader.close()
        }        
    }
    
    void stopActor(String name, Actor actor) {
        log.info "Stopping $name"
        actor << "stop"
        actor.join()
    }
}
