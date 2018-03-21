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
    Map<String,CompactReadPair> buffer 
    
    Actor consumer
    
    int received = 0
    
    int paired
    
    int chimeric
    
    Regions regions
    
    String debugRead = null // "SN7001291:342:HFMC7BCXX:2:1102:7840:64817"
    
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
        
        final String readName = record.readName
        CompactReadPair pair = buffer[readName]
        if(pair) {
            consumer << [pair, record]
            buffer.remove(readName)
            paired += 2
            return
        }
        
        // Pair not found - create a new one
        pair = new CompactReadPair(record)
        
        // Is it or its mate located in the desired regions?
        if(regions != null) 
            if(readName == debugRead) {
                log.info "Debug read: $record"
            }
            
            if(notInRegion(pair, readName == debugRead)) {
                if(readName == debugRead)
                    log.info "$readName filtered out"
                return
            }
        
        if(pair.chimeric)
            ++this.chimeric
            
        buffer[readName] = pair
    }
    
    boolean notInRegion(CompactReadPair pair, boolean debug) {
        
        if(pair.unmapped) {
            return false
        }
        
        final String chr1 = pair.r1ReferenceName
        final int r1p = pair.r1AlignmentStart
        final int len = pair.readLength
        if(regions.overlaps(chr1, r1p , r1p+len))
            return false
        
        final String chr2 = pair.r2ReferenceName?:chr1
        final int r2p = pair.r2AlignmentStart
        if(regions.overlaps(chr2, r2p , r2p+len)) {
            return false 
        }
        
        if(debug)
            log.info "Regions do not overlap $chr2:$r2p-${r2p+len}"
       
        return true
    }
    
    String toString() {
        [received,paired,buffer.size()].collect { Utils.human(it) }.join(',')
    }
}