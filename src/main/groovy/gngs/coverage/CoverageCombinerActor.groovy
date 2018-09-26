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
package gngs.coverage

import java.util.concurrent.atomic.AtomicInteger

import gngs.*
import groovy.transform.CompileStatic
import groovy.util.logging.Log
import groovyx.gpars.actor.Actor
import groovyx.gpars.actor.DefaultActor
import htsjdk.samtools.SAMRecordIterator

@Log
class CoverageCombinerActor extends RegulatingActor<SampleReadCount> {
    
    /**
     * How many samples to expect
     */
    int numSamples = 0
    
    String chr
    
    long pos = -1
    
    TreeMap<Long,Map<String,Integer>> counts = new TreeMap()
    
    AtomicInteger countSize = new AtomicInteger()
    
    AtomicInteger pending = new AtomicInteger()
    
    AtomicInteger processed = new AtomicInteger()
    
    Map<String,AtomicInteger> pendingCounts = Collections.synchronizedMap([:])
    
    ProgressCounter progress = new ProgressCounter(withRate: true, timeInterval: 1000, lineInterval: 200, log:log, 
        extra: { "Combining $chr:$pos (${XPos.parsePos(pos).startString()}) with ${counts[pos]?.size()?:0}/$numSamples reported, count buffer size=${countSize} (Samples = ${counts[pos]})" })
    
    
    CoverageCombinerActor(RegulatingActor printer, int numSamples) {
        super(printer, 20000, 200000)
        this.progress = progress
        this.numSamples = numSamples
    }
    
    void add(SampleReadCount sampleCount) {
        pending.incrementAndGet()
        this << sampleCount
    }
    
    @CompileStatic
    void process(SampleReadCount count) {
        
        pending.decrementAndGet()
        
        this.chr = count.chr
        
        Long xpos = XPos.computePos(count.chr, count.pos)
        if(pos == -1)
            pos = xpos
        
        Map<String,Integer> positionCounts = this.counts[xpos]
        if(positionCounts == null) {
            this.counts[xpos] = positionCounts = new HashMap()
            countSize.incrementAndGet()
        }
        
        positionCounts[count.sample] = count.reads
        progress.count()
        
        pendingCounts.get(count.sample, new AtomicInteger(0)).getAndIncrement()
        
        if(pos == xpos) {
            if(positionCounts.size() == numSamples) {
                        
                counts.pollFirstEntry()
                
                countSize.decrementAndGet()
                
                if(!counts.isEmpty())
                    pos = counts.firstKey()
                else
                    pos = -1
                    
                sendDownstream([region: count.target, chr: count.chr, pos: count.pos, counts: positionCounts])
           }
        }
    }
     
    @CompileStatic
    void processBAM(SAM bam, Regions scanRegions) {
        CoverageCalculatorActor.processBAM(bam, scanRegions, this)
    }
    
    String toString() {
        "CoverageCombinerActor(${XPos.parsePos(pos).startString()})"
    }
}

