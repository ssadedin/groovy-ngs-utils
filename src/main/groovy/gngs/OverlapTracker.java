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
package gngs;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.TreeSet;

import gngs.coverage.ReadRange;
import htsjdk.samtools.SAMRecord;

/**
 * Simple class to track ranges and remove them when they no longer overlap a given 
 * position.
 * <p>
 * This turns out to be a hotspot for performance in coverage depth 
 * calculations.
 * 
 * @author Simon Sadedin
 */
public final class OverlapTracker {
    
    public final List<int[]> reads = new LinkedList<int[]>();
    
    public final void add(SAMRecord record) {
        reads.add(new int[] { record.getAlignmentStart(), record.getAlignmentEnd()});
    }
    
    public final void add(ReadRange r) {
        reads.add(new int [] { r.alignmentStart, r.alignmentEnd });
    }
    
    public int size() {
        return this.reads.size();
    }
    
    public void clear() {
        this.reads.clear();
    }
    
    public void removeNonOverlaps(final int pos) {
        final Iterator<int[]> iter = reads.iterator();
        while(iter.hasNext()) {
            int[] readRange = iter.next();
            if(pos>readRange[1])
                iter.remove();
        }
    }
}
