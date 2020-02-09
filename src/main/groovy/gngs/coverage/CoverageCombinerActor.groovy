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
import groovy.transform.ToString
import groovy.util.logging.Log
import groovyx.gpars.actor.Actor
import groovyx.gpars.actor.DefaultActor
import htsjdk.samtools.SAMRecordIterator
import htsjdk.tribble.readers.TabixReader

@CompileStatic
@ToString
class PositionCounts {
    final Region region
    final String chr
    final int pos
    final Map<String,Integer> counts

    public PositionCounts(final Region region, final int pos, final Map<String, Integer> counts) {
        super();
        this.region = region;
        this.chr = region.chr;
        this.pos = pos;
        this.counts = counts;
    }    
}

@Log
class CoverageCombinerActor extends RegulatingActor<SampleReadCount> {
    
    /**
     * How many samples to expect
     */
    int numSamples = 0
    
    String chr
    
    long pos = -1

    TreeMap<Long,Map<String,Integer>> counts = new TreeMap()
    
    final AtomicInteger countSize = new AtomicInteger()
    
    final AtomicInteger pending = new AtomicInteger()
    
    final AtomicInteger processed = new AtomicInteger()

    boolean countFragments = true
    
    final Map<String,AtomicInteger> pendingCounts = Collections.synchronizedMap([:])
    
    ProgressCounter progress = new ProgressCounter(withRate: true, timeInterval: 1000, lineInterval: 200, log:log, 
        extra: this.&statusMessage)
    
    Map<String,Long> samplePositions = new HashMap()
    
    String statusMessage() {
        "Combining $chr:$pos (${XPos.parsePos(pos).startString()}) with ${counts[pos]?.size()?:0}/$numSamples reported, count buffer size=${countSize} (Samples = ${counts[pos]})"  
    }
    
    CoverageCombinerActor(RegulatingActor printer, int numSamples, boolean countFragments = true) {
        super(printer, 20000, 200000)
        this.progress = progress
        this.numSamples = numSamples
        this.countFragments = countFragments
    }
    
    void add(SampleReadCount sampleCount) {
        pending.incrementAndGet()
        this << sampleCount
    }
    
    @CompileStatic
    void process(final SampleReadCount count) {
        
        final String sample = count.sample
        
        pending.decrementAndGet()
        
        this.chr = count.chr
        
        final Long xpos = XPos.computePos(count.chr, count.pos)
        if(pos == -1)
            pos = xpos
        
        final Map<String,Integer> positionCounts = getOrInitCounts(xpos)
        positionCounts.put(sample,count.reads)
        
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
                    
                PositionCounts posCounts = new PositionCounts(count.target, count.pos, positionCounts)
//                log.info "Send: $posCounts"
                sendDownstream(posCounts)
           }
        } 
    }
    
    @Override
    void onEnd() {
        progress.end()
        log.info "CoverageCombiner stopped at $pos"
    }

    @CompileStatic
    private Map getOrInitCounts(long xpos) {
        Map<String,Integer> positionCounts = this.counts[xpos]
        if(positionCounts == null) {
            this.counts[xpos] = positionCounts = new HashMap()
            countSize.incrementAndGet()
        }
        return positionCounts
    }
     
    @CompileStatic
    void processBAM(SAM bam, Regions scanRegions, int minMQ, int downsampleWindow=0, int subsample=1, String sample = null, boolean countFragments = true) {
        
        RegulatingActor<SampleReadCount> coverageSink = this
        boolean stopSink = false
        if(downsampleWindow>1) {
            log.info "Adding coverage downsampler because downsampleWindow=$downsampleWindow, subsample=$subsample"
            coverageSink = new CoverageDownsampler(this, downsampleWindow, subsample)
            coverageSink.start()
            stopSink = true
        }
       
        CoverageCalculatorActor.processBAM(bam, scanRegions, coverageSink, minMQ, sample, countFragments)

        if(stopSink) {
            coverageSink.sendStop()
            coverageSink.join()
        }
    }
    
    @CompileStatic
    void processTabix(String sample, TabixReader reader, Regions scanRegions, int downsampleWindow=0, int subsampleBy) {
        
        RegulatingActor<SampleReadCount> coverageSink = this
        boolean stopSink = false
        if(downsampleWindow>1) {
            coverageSink = new CoverageDownsampler(this, downsampleWindow, subsampleBy)
            coverageSink.start()
            stopSink = true
        }
       
        AtomicInteger downstreamCount = new AtomicInteger(0)
        
        List<String> chrs = scanRegions*.chr.unique()
        int count = 0
        BatchedAcknowledgeableMessage bam = new BatchedAcknowledgeableMessage(downstreamCount, 10)
        for(String chr in chrs) {
            int start = Math.max(0, scanRegions.index[chr].ranges.firstKey() - 1000)
            int end = scanRegions.index[chr].ranges.lastKey() + 1000
            log.info "Scan $chr from $start to $end (overlapping ${scanRegions.numberOfRanges} regions)"
            Region scanRegion = new Region(chr, start, end)
            TabixReader.Iterator iter = reader.query(chr, start-2, end+2)
            String line
            
            Regions overlapRegions = scanRegions.collect { Region r -> r.widen(-1,0) } as Regions
            boolean emitting = false
            int lastPos = -1
            int debugPos = 237551483
            while((line = iter.next()) != null) {
                List<String> fields = line.tokenize('\t')
                
                final int pos = Integer.parseInt(fields[1])
                final int readCount = Integer.parseInt(fields[2])
                
                if(overlapRegions.overlaps(chr, pos, pos)) {
                    SampleReadCount countInfo = new SampleReadCount(
                        scanRegion,
                        chr,
                        pos,
                        readCount,
                        sample
                    )
                    bam.batchTo(countInfo, coverageSink)
                    emitting = true
                }
                ++count
            }
            bam.flush(this)
            if(stopSink) {
                coverageSink.sendStop()
                coverageSink.join()
            }
        }
    }
    
    String toString() {
        "CoverageCombinerActor(${XPos.parsePos(pos).startString()})"
    }
}

