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
 * Discovers blocks of low coverage in a file output in BEDTools Coverage format,
 * or alternatively {@link MultiCov} output.
 * <p>
 * Each low coverage block is output as a {@link CoverageBlock} object which 
 * has statistics about the coverage within the block. As {@link CoverageBlock} extends
 * {@link Region}, you can easily turn this into a set of Regions to perform more
 * calculations or statistics.
 * <p>
 * To use this class, create the object first by passing a file to compute coverage for,
 * then call the {@link #calculate} method to actually calculate the coverage blocks:
 * <pre>
 * CoverageGaps gaps = new CoverageGaps("test.cov.gz")
 * gaps.calculate()
 * println("There are " + gaps.blocks.size() + " low coverage regions")
 * </pre>
 * @author Simon Sadedin
 */
@Log
class CoverageGaps {
    
    /**
     * Use a large buffer for very efficient reading
     */
    int BUFFER_SIZE_BYTES = 1024*1024*20
    
    /**
     * Only include gaps that are at least this size
     */
    int minRegionWidth = 0
    
    String coverageFilePath
    
    int threshold = 20
    
    int lineCount = 0
    int blockCount = 0
    int totalBP = 0
    
    
    SummaryStatistics coverageStats = new SummaryStatistics()
    CoverageStats coveragePercentiles = new CoverageStats(1000)
    
    /**
     * Discovered gaps (the results) are stored here for access after calling {@link #calculate}
     */
    List<CoverageBlock> blocks = []
    
    CoverageBlock block = null
    
    Regions regions = null
    
    int offset = 0
    
    RegulatingActor<CoverageBlock> gapProcessor = null
    
    CoverageGaps(String coverageFilePath) {
        this.coverageFilePath = coverageFilePath
    }
    
    void calculate() {
		Utils.reader(coverageFilePath) { 
            calculateFromBEDTools(it)
		}
    }
    
    void calculateMultiCov() {
		Utils.reader(coverageFilePath) { 
            calculateFromMultiCov(it)
        }
    }
    
    @CompileStatic
    void calculateFromMultiCov(Reader reader) {
        String chr = null
        int start = -1
        int end = -1
        
        String line
        int lastPos = -1
        int covField = -1
        boolean first = true
        ProgressCounter progress = createProgress()
        while((line = reader.readLine()) != null) {
            progress.count()
            List<String> fields = line.tokenize('\t')
            chr = fields[0]
            
            // Use the first line to check if the coverage is in the 3rd field or the 4th
            if(first) {
                first = false
                covField = fields.size() == 3 ? 2 : 3
            }
            
            int pos = Integer.parseInt(fields[1])
            int cov = Integer.parseInt(fields[covField])
            if(pos != lastPos + 1) {
                start = pos
            }
            processLine(chr, start, pos, cov, null)
            lastPos = pos
        }
        progress.end()
    }
    
    @CompileStatic
    void calculateFromBEDTools(Reader reader) {
        
        ProgressCounter progress = createProgress()
            
        int pos = 0
        
        reader.eachLine { String line ->
            
            progress.count()
            
            List<String> fields = line.tokenize('\t')
            
            int cov = fields[-1].toInteger()
            String chr = fields[0]
            int start = fields[1].toInteger()
            int end = fields[2].toInteger()
            int offset = fields[-2].toInteger()
            
            String id = null
            if(fields.size()>5)
                id = fields[-3] 
            
            processLine(chr, start, start+offset, cov, id)
        }
        
        // we might still have a block
        if(block && (block.end - block.start + 1 >= minRegionWidth)) {
            log.info("Final block processed")
            outputBlock(pos)
            block = null;
        }
        
        progress.end()
    }

    private ProgressCounter createProgress() {
        ProgressCounter progress =
                new ProgressCounter(withTime: true, withRate:true, timeInterval:10000, log:log, extra: {
                    "Region $region, ${blocks.size()} gaps after scanning ${totalBP} bp"
                })
        return progress
    }
    
    String region
    
    @CompileStatic
    void processLine(String chr, int start, int pos, int cov, String id) {
        coverageStats.addValue(cov)
        coveragePercentiles.addValue(cov)
        region = "$chr:$start"
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
                if(id != null)
                    block.id = id
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
    
    @CompileStatic
    void outputBlock(int endPos=-1) {
        blockCount++
        if(endPos>=0)
            block.end = endPos
        
        if(block.end < block.start)
            assert false
            
        blocks.add(block)
        
        if(this.gapProcessor != null)
            this.gapProcessor.sendTo(block)
            
        block = null
    }
    
    @CompileStatic
    Object asType(Class clazz) {
        if(clazz == Regions) {
            new Regions((Iterable<IRegion>)blocks.collect { it as Region })
        }
    }
}


