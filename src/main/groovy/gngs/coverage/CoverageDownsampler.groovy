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

import gngs.RegulatingActor
import graxxia.RollingStats
import groovy.transform.CompileStatic
import groovy.util.logging.Log

/**
 * Accepts SampleReadCount objects as messages and calculates a moving window
 * of their values, emitting the window mean instead of the original
 * values.
 * <p>
 * In addition to calculating the moving window, this class can also subsample
 * the results so that only every 1 out of N values is emitted. The size of the
 * moving window and the subsampling factor can be controlled independently.
 * 
 * @author simon.sadedin
 */
@CompileStatic
@Log
class CoverageDownsampler extends RegulatingActor<SampleReadCount> {
    
    RollingStats stats
    
    final int downsampleBy
    
    final int subsampleBy
    
    private int index = 0
    
    private int countEmitted = 0
    
    private int lastPos = -1
    
    private SampleReadCount lastCounts
    
    CoverageDownsampler(RegulatingActor<SampleReadCount> downstream, int windowSize, int subsampleBy) {
        super(downstream, 5000,10000)
        this.downsampleBy = windowSize
        this.subsampleBy = subsampleBy
    }

    @Override
    public void process(final SampleReadCount counts) {
        
        if(lastPos != counts.pos-1) {
            newRegion()
        }
        
        lastPos = counts.pos
        lastCounts = counts 
        
        stats.addValue(counts.reads)
        
        ++index
        if(index % subsampleBy == 0) {
            emit(counts)
            ++countEmitted
        }
    }

    private void emit(SampleReadCount counts) {
        int meanEstimate = (int)Math.round(stats.stableMean)
        SampleReadCount downsampled =
                        new SampleReadCount(counts.target, counts.chr, counts.pos, meanEstimate, counts.sample)
        sendDownstream(downsampled)
    }
    
    private void newRegion() {
        if(countEmitted == 0) {
            if(lastCounts != null) {
                log.info "Emitting last counts for ${lastCounts.target} because no counts emitted from downsampler"
                emit(lastCounts)
            }
        }            
        this.countEmitted = 0
        this.lastCounts = null
        this.stats = new RollingStats(downsampleBy)
        this.index = 0
    }
}
