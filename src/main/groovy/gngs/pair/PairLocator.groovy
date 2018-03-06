package gngs.pair

import java.io.Writer

import gngs.*

import groovy.transform.CompileStatic
import groovy.util.logging.Log
import groovyx.gpars.actor.Actor
import groovyx.gpars.actor.DefaultActor
import htsjdk.samtools.SAMRecord

@CompileStatic
@Log
class PairLocator extends DefaultActor {
    
    /**
     * Read pairs we are indexing, keyed on name
     */
    Map<String,SAMRecordPair> buffer 
    
    Actor consumer
    
    int received = 0
    
    int paired
    
    int chimeric
    
    Regions regions
    
    PairLocator(Actor consumer) {
        this.consumer = consumer
        this.buffer = new HashMap(200000)
    }
    
    void act() {
        loop {
            react { msg ->
                assert msg != null
                
                if(msg == "stop") {
                    terminate()
                }
                else {
                    processRead((SAMRecord)msg)
                }
            }
        }
    }
    
    void processRead(SAMRecord record) {
        
        ++received
        
        if(record.secondaryOrSupplementary)
            return
        
        SAMRecordPair pair = buffer[record.readName]
        if(pair) {
            pair.r2 = record
            consumer << pair
            buffer.remove(record.readName)
            paired += 2
            return
        }
        
        // Pair not found - create a new one
        pair = new SAMRecordPair(r1: record)
        
        // Is it or its mate located in the desired regions?
        if(regions != null) 
            if(notInRegion(pair))
                return
        
        if(pair.chimeric)
            ++this.chimeric
            
        buffer[record.readName] = pair
    }
    
    boolean notInRegion(SAMRecordPair pair) {
        
        if(pair.unmapped) {
            return false
        }
        
        String chr1 = pair.r1ReferenceName
        int r1p = pair.r1Pos
        int len = pair.readLength 
        if(regions.overlaps(chr1, r1p , r1p+len))
            return false
        
        String chr2 = pair.r2ReferenceName
        int r2p = pair.r2Pos
        if(regions.overlaps(chr2, r2p , r2p+len)) {
            return false 
        }
       
        return true
    }
    
    String toString() {
        [received,paired,buffer.size()].collect { Utils.human(it) }.join(',')
    }
}