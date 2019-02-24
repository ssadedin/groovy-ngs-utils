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

@CompileStatic
class CoverageDownsampler extends RegulatingActor<SampleReadCount> {
    
    RollingStats stats
    
    final int downsampleBy
    
    int index = 0
    
    int lastPos = -1
    
    CoverageDownsampler(RegulatingActor<SampleReadCount> downstream, int windowSize) {
        super(downstream, 5000,10000)
        this.downsampleBy = windowSize
    }

    @Override
    public void process(final SampleReadCount counts) {
        
        if(lastPos != counts.pos-1)
            this.stats = new RollingStats(downsampleBy)
        
        lastPos = counts.pos
        
        stats.addValue(counts.reads)
        
        SampleReadCount downsampled = new SampleReadCount(counts.target, counts.chr, counts.pos, (int)Math.round(stats.stableMean), counts.sample) 
        sendDownstream(downsampled)
    }
}
