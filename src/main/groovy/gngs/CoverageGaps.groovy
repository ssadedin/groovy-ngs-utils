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

import groovy.lang.IntRange
import groovy.transform.CompileStatic
import groovy.util.logging.Log

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics
import org.apache.commons.math3.stat.descriptive.SummaryStatistics

/**
 * A run of contiguous bases within a band of coverage depth
 * (eg: a run of low coverage bases).
 * <p>
 * Note that the start and end values are <b>inclusive</b>. That is,
 * each base including both start and end values is below coverage 
 * threshold.
 * 
 * @author Simon Sadedin
 */
@CompileStatic
class CoverageBlock implements IRegion {
    String region
    String chr
    int start
    int end
    String id
    DescriptiveStatistics stats = new DescriptiveStatistics()
    
    Object asType(Class clazz) {
        if(clazz == Region) {
            return new Region(chr, new GRange(start,end, this))
        }
    }

    @Override
    public IntRange getRange() {
        return start..end
    }
    
    int size() {
        return end - start + 1
    }
}

/**
 * Discovers blocks of low coverage in a file output in BEDTools Coverage format.
 * <p>
 * Each low coverage block is output as a {@link CoverageBlock} object which 
 * has statistics about the coverage within the block. 
 * <p>
 * To use this class, create the object first by passing a file to compute coverage for,
 * then call the {@link #calculate} method to actually calculate the coverage blocks:
 * <pre>
 * CoverageGaps gaps = new CoverageGaps("test.cov.gz")
 * gaps.calculate()
 * </pre>
 * @author Simon Sadedin
 */
@Log
class CoverageGaps {
    
    /**
     * Use a large buffer for very efficient reading
     */
    int BUFFER_SIZE_BYTES = 1024*1024*20
    
    int minRegionWidth = 0
    
    String coverageFilePath
    
    int threshold = 20
    
    int lineCount = 0
    int blockCount = 0
    int totalBP = 0
    
    
    SummaryStatistics coverageStats = new SummaryStatistics()
    CoverageStats coveragePercentiles = new CoverageStats(1000)
    
    List<CoverageBlock> blocks = []
    
    CoverageBlock block = null
    
    int offset = 0
    
    CoverageGaps(String coverageFilePath) {
        this.coverageFilePath = coverageFilePath
    }
    
    void calculate() {
		InputStream stream = new BufferedInputStream(new FileInputStream(coverageFilePath), BUFFER_SIZE_BYTES)
		if(coverageFilePath.endsWith(".gz")) 
			stream = new java.util.zip.GZIPInputStream(stream) 

        calculateFromStream(stream)
    }
    
    @CompileStatic
    void calculateFromStream(InputStream stream) {
        
        ProgressCounter progress = 
            new ProgressCounter(withTime: true, withRate:true, timeInterval:10000, log:log, extra: {
            "${blocks.size()} gaps after scanning ${totalBP} bp"
        })
        
        int pos = 0
        
        stream.eachLine { String line ->
            
            progress.count()
            
            List<String> fields = line.tokenize('\t')
            
            int cov = fields[-1].toInteger()
            String chr = fields[0]
            int start = fields[1].toInteger()
            int end = fields[2].toInteger()
            int offset = fields[-2].toInteger()
            
            coverageStats.addValue(cov)
            coveragePercentiles.addValue(cov)
            pos = start + offset
            String region = "$chr:$start"
            ++totalBP

            if(block && block.region != region) { // end of region
                if (block.end - block.start + 1 >= minRegionWidth) { // only write if long enough
                    outputBlock()
                }
                else {
                    block = null; // forget this (too short) block
                }
            }

            if(cov < threshold) {
                if(!block)  {
                   block = new CoverageBlock(chr:chr, region:region, start:pos)
                   if(fields.size()>5)
                       block.id = fields[-3]
                }
                block.stats.addValue(cov)
                block.end = pos // note that end is the last position that is a gap
            }
            else { // coverage is ok
                if(block && (block.end - block.start + 1 >= minRegionWidth)) {
                    outputBlock()
                }
                else { // forget any block as it's too short
                    block = null;
                }
            }
        }
        
        // we might still have a block
        if(block && (block.end - block.start + 1 >= minRegionWidth)) {
            log.info("Final block processed")
            outputBlock(pos)
            block = null;
        }
        
        progress.end()
    }
    
    @CompileStatic
    void outputBlock(int endPos=-1) {
        blockCount++
        if(endPos>=0)
            block.end = endPos
        
        if(block.end < block.start)
            assert false
            
        blocks.add(block)
        block = null
    }
    
    @CompileStatic
    Object asType(Class clazz) {
        if(clazz == Regions) {
            new Regions((Iterable<IRegion>)blocks.collect { it as Region })
        }
    }
}


