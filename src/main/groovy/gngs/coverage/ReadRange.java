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
package gngs.coverage;

import htsjdk.samtools.SAMRecord;

public class ReadRange {
    
    public ReadRange(SAMRecord r) {
        this(r, true); // by default count fragments, not reads
    }
    
    public ReadRange(SAMRecord r, boolean countFragments) {
        
        this.referenceIndex = r.getReferenceIndex();
        this.alignmentStart = r.getAlignmentStart();
        int alignmentEnd = r.getAlignmentEnd();
        final int mateStart = r.getMateAlignmentStart();
        
        // Avoid double counting if two reads from the same fragment will overlap
        if(countFragments && (r.getFirstOfPairFlag() && mateStart >= this.alignmentStart && mateStart <= alignmentEnd)) {
            alignmentEnd = mateStart;
        }
        this.referenceName = r.getReferenceName();
        this.alignmentEnd = alignmentEnd;
    }
    
    final public int referenceIndex;
    
    final public int alignmentStart;
    
    final public int alignmentEnd;
    
    final public String referenceName;
    
    public String toString() {
        return String.format("%s:%d-%d", referenceName, alignmentStart, alignmentEnd);
    }
}
