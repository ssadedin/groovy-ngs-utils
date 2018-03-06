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
    
    PairLocator chimericLocator
    
    PairFormatter formatter 
    
    Regions regions
    
    List<Actor> actors = []
    
    int chimeric
    
    int shardId = -1
    
    int shardSize = 0
    
    int numLocators 
    
    PairScanner(Writer writer, int numLocators, Regions regions = null) {
        this.pairWriter = new PairWriter(writer)
        this.formatter = new PairFormatter(1000 * 1000, pairWriter)
        this.chimericLocator = new PairLocator(formatter)
        this.regions = regions
        this.numLocators = numLocators
        progress.log = log
    }
    
    void initLocators() {
        int locatorsCreated = 0
        this.locators = []
        
        // Note: since we may be running in sharded mode, certain locator positions
        // may never get used. We counter this by checking at each position
        // if the locator will be used and we keep creating locators until we 
        // have the requisite number at each of the shard positions in our array
        for(int i=0; locatorsCreated < numLocators; ++i) {
            if(shardId<0 || ((i%shardSize) == shardId)) {
                ++locatorsCreated
                PairLocator pl = new PairLocator(formatter)
                if(this.regions)
                    pl.regions = new Regions(this.regions)        
                this.locators << pl
                this.locatorIndex << pl
            }
            else {
                this.locatorIndex.add(null)
            }
        }
        
        this.actors = locators + [
            chimericLocator,
            formatter,
            pairWriter,
        ]
        this.actors*.start()
    }
    
    @CompileStatic
    void scan(SAM bam) {
        log.info "Beginning scan of $bam.samFile"
        running = this
        this.initLocators()
        try {
            this.scanBAM(bam)
        }
        finally {
            log.info "Stopping parallel threads ..."
            
            locators.eachWithIndex { a, i -> stopActor("Locator $i", a) }
            
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
                    if(locator != null)
                        locator << read 
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
