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
import java.util.logging.Logger

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics
import org.apache.commons.math3.stat.descriptive.SummaryStatistics
import groovy.transform.CompileStatic
import groovy.util.logging.Log
import groovyx.gpars.actor.Actor
import groovyx.gpars.actor.DefaultActor
import htsjdk.samtools.SamReader
import htsjdk.samtools.SAMFormatException
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SAMRecordIterator

@CompileStatic
class ResolveRequest {
    
    long scheduleTimeMs = System.currentTimeMillis()
    
    SAMRecordPair pair
}

/**
 * An actor that can accept {@link ResolveRequest} messages to asynchronously 
 * resolve the missing mates of read pairs.
 * <p>
 * <b>Note</b>: this class i superceded by the {@link AsyncMateResolveWorker} class,
 * and eventually should be refactored into that class. This is because the 
 * a "pull" model where workers pull work requests from a queue seems to be a better
 * fit than a push model where the sender decides the priority.
 * 
 * @author Simon Sadedin
 */
class ReadMateLookupActor extends DefaultActor implements Closeable {
    
    private static Logger log = Logger.getLogger("ReadMateLookupActor")
    
    SamReader reader
    
    int resolved = 0
    
    Closure resolvedCallback
    
    /**
     * Mate positions in buffer
     */
    MissingMateIndex mateIndex
    
    /**
     * Count of how many lookaheads we have done
     */
    int lookAheads = 0
    
    int total = 0
    
    Map<String,SummaryStatistics> stats = [
        time: new SummaryStatistics(),
        latency: new SummaryStatistics(),
        scanSize: new SummaryStatistics(),
        efficiency:new SummaryStatistics() 
    ]
    
    /**
     * How many unresolved pairs got resolved by the ahead lookups we have done
     */
    int lookAheadResolutions = 0
    
    int presolved = 0
    
    public ReadMateLookupActor(SamReader reader, MissingMateIndex mateIndex, Closure resolvedCallback) {
        super();
        this.reader = reader;
        this.resolvedCallback = resolvedCallback
        this.mateIndex = mateIndex
    }

    @Override
    void act() {
        loop { 
            react { msg ->
                if(msg instanceof ResolveRequest) {
                    assert msg != null
                    resolvePair(msg)
                }
                else
                if(msg == "stop") {
                    terminate()
                }
            }
        }
    }
    
    String debugRead = "H233NCCXX150123:6:2216:5912:67112"
    
    @CompileStatic
    void resolvePair(ResolveRequest request) {
        SAMRecordPair pair = request.pair
        
        stats.latency.addValue(System.currentTimeMillis() - request.scheduleTimeMs)
        
        // It could be that the mate index was cleared of this pair while we
        // were waiting to process it. If so, we should skip it
        assert mateIndex != null
        assert pair != null
        assert resolvedCallback != null
        
        ++total
        
        SAMRecordPair check = mateIndex.getPair(pair.readName)
        if(check == null) {
            ++presolved
            return
        }
        
        if(pair.hasBothReads()) {
            ++presolved
            return
        }
        
        ScannableRegion scannableRegion = mateIndex.searchForScannableRegion(pair, 500, 50)
        
        Region regionToScan = scannableRegion.region
        
        if(regionToScan == null)
            regionToScan = new Region(pair.nonMissingRead.mateReferenceName, pair.missingReadPosition-50, pair.missingReadPosition+50)
        
        resolvePair(pair, regionToScan)
    }
    
     @CompileStatic
    void resolvePair(SAMRecordPair pair, Region regionToScan) {
        
        if(pair.r1 != null) {
            if(regionToScan != null) {
//                if(debugRead != null && pair.readName == debugRead) {
//                    log.info "Scanning $pair.readName with region scan"
//                }
                
                SAMRecord foundMate = queryMateByRegionScan(pair.r1, regionToScan)
                if(foundMate != null) {
                    pair.r2 = foundMate
                }
                else
                    log.info "Failed to find mate for $pair.readName via region scan of $regionToScan even though mate pos = ${pair.r1.mateReferenceName}:${pair.r1.mateAlignmentStart}"
            }
            else {
                pair.r2 = queryMate(pair.r1)
            }
            
            if(pair.r2 != null) {
                ++resolved
            }
            resolvedCallback.call(pair)
        }
        else 
        if(pair.r2 != null) {
            if(regionToScan != null) {
                SAMRecord foundMate = queryMateByRegionScan(pair.r2, regionToScan)
                if(foundMate != null) {
                    pair.r1 = foundMate
                }
            }
            else {
                pair.r1 = queryMate(pair.r2)
            }
            
            if(pair.r1 != null) {
                ++resolved
            }
            resolvedCallback.call(pair)
        }            
    }
    
    @CompileStatic
    SAMRecord queryMate(SAMRecord r1) {
        
        
        try {
            return reader.queryMate(r1)
        }
        catch(SAMFormatException sfe) {
            // ignore
            // happens due to BWA secondary alignments
        }
        return null
    }
    
    @CompileStatic
    SAMRecord queryMateByRegionScan(SAMRecord record, Region region) {
        
        long scanStartTimeMs = System.currentTimeMillis()
        
        ++lookAheads
        
        SAMRecordIterator i = this.reader.query(region.chr, 
                                                       region.from, 
                                                       region.to, false)
        SAMRecord mate = null
        int readCount = 0
        int resolvedCount = 0
        try {
            while(i.hasNext()) {
                ++readCount
                SAMRecord potentialMate = i.next()
                if(!potentialMate)
                    continue
                    
                if(potentialMate.readName == record.readName) {
                    mate = potentialMate
                    ++resolvedCount
                }
                else {
                    SAMRecordPair pair = mateIndex.getPair(potentialMate.readName)
                    if(pair) { 
                        if(pair.r2 == null) {
                            pair.r2 = potentialMate
                            mateIndex.resolve(potentialMate.alignmentStart, pair)
                            ++lookAheadResolutions
                        }
                    }
                    else { 
                        // We scanned a read that was not being looked for (yet)
                        // However it could be that we *will* look for it
                        // If the read is in the same chromosome, chances are, finding it will be
                        // done by the (efficient) linear scan. But if the read is
                        // chimeric we will probably be forced to do an inefficient random lookup
                        // Therefore we will choose to index chimeric reads in case they show up 
                        // later as useful
                        if(potentialMate.referenceIndex != potentialMate.mateReferenceIndex && potentialMate.mateReferenceIndex == record.referenceIndex) {
                            mateIndex.add(new SAMRecordPair(r1:potentialMate))
                        }
                    }
                }
            }
            return mate
        }
        finally {
            addStats(System.currentTimeMillis()-scanStartTimeMs, readCount, resolvedCount, region)
            i.close()
        }
    }
    
    public void addStats(long timeMs, int readCount, int resolvedCount, Region region) {
        
        stats.time.addValue(timeMs)
        stats.scanSize.addValue(region.size())
        stats.efficiency.addValue(resolvedCount/(readCount+1))
        
        if(lookAheads % 500 == 0) {
            log.info "${this.hashCode()} resolved ${resolved}/ ${total}, presolved ${presolved}, latency=${Math.round(stats.latency.mean)}, lookaheads $lookAheads (mean time=${Math.round(stats.time.mean)}ms), " +
                     "scanSize=${Utils.humanBp(Math.round(stats.scanSize.mean))} " +
                     "efficiency=${Utils.perc(stats.efficiency.mean)}"
        }
    }

    @Override
    public void close() {
        Utils.closeQuietly(reader) 
    }
}
