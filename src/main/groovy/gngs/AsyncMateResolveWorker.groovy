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

import groovy.transform.CompileStatic
import groovy.util.logging.Log

/**
 * A background task that can read from a queue of unpaired reads and resolve their
 * pairs in the background, as well as resolve other reads that are in the same
 * region as the pair that is resolved.
 * <p>
 * This class interoperates closely with the {@link MissingMateIndex} class, which
 * maintains a priority queue (see {@link MissingMateIndex#resolveQueue}) of reads that
 * are missing mates which need resolving. 
 * <p>
 * **Note**: due to the blocked nature of BAM format files, it is highly inefficient to 
 * read just one read from a block, but comparatively extremely efficient to read all of them.
 * So the logic applied is that this class takes a pair to resolve from the queue, but then
 * it doesn't just resolve that pair, it will resolve any pair that it happens to read the 
 * mate for that also needs resolution and which is in the same region.
 * 
 * @author Simon Sadedin
 */
@Log
@CompileStatic
class AsyncMateResolveWorker implements Runnable {
    
    MissingMateIndex index
    
    ReadMateLookupActor lookupActor
    
    ProgressCounter progress = new ProgressCounter(withTime:true, extra: {
        "unresolvable=$unresolvable"
    })
    
    int unresolvable = 0
    
    /**
     * Main loop - runs forever
     */
    void run() {
        while(true) {
            
            RegionScan scan
            SAMRecordPair pair
            synchronized(index) {
                if(index.resolveQueue.isEmpty())
                    index.wait()
                 else {
                     pair = index.resolveQueue.poll()
                 }
                 
                 if(!needsResolution(pair))
                     continue
                     
                 scan = index.allocateScannableRegion(pair, 500, 1)
                 if(scan == null) {
                     // we got a pair back but its mate was not in the mate index
                     // the only mates not indexed are non-major chromosomes, etc
                     ++unresolvable
                     continue
                 }
             }
             
             progress.count()
             
             Region region = scan.region.region
             lookupActor.resolvePair(pair, region)
             
             scan.markUnpaired()
        }
    }
    
    /**
     * Check if a read pair is still in a state where the mate needs to be resolved
     * 
     * @param pair  the read pair to check
     * @return      true if this read pair needs resolution still
     */
    boolean needsResolution(SAMRecordPair pair) {
         if(pair == null)
             return false
                     
         if(index.getPair(pair.readName) == null) // must have been resolved
             return false
                 
         if(pair.hasBothReads())
             return false
                     
         // another async resolver could have decided this is unpairable
         if(pair.flags.unpaired) 
             return false
                     
         return true
    }
}