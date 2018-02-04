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
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SAMTagUtil
 
/**
 * Models a read pair, where one of the reads might be missing.
 * <p>
 * In a set of paired reads there are a number of attributes that are recorded
 * in both pairs. So it is convenient to be able to write code in an abstract
 * way that can query these attributes regardless of which read of the pair
 * (or both) are actually available.
 * 
 * @author Simon Sadedin
 */
class SAMRecordPair implements Comparable {
    
    SAMRecord r1
    
    SAMRecord r2
    
    Map<String,Object> flags = Collections.synchronizedMap([:])
    
    @CompileStatic
    boolean hasBothReads() {
        r2 != null && r1 != null
    }
    
    @CompileStatic
    boolean isChimeric() {
        if(r1 != null)
            return r1.referenceIndex != r1.mateReferenceIndex
        else
        if(r2 != null)
            return r2.referenceIndex != r2.mateReferenceIndex
        else
            return false // both reads null?!
    }
    
    @CompileStatic
    void setTag(String name, String value) {
        r1.setAttribute(name, value)
        r2.setAttribute(name, value)
    }
    
    static short RG_TAG = SAMTagUtil.getSingleton().RG
    
    /**
     * Do not static compile this as it relies on access to protected function
     *
     * @param rgId
     */
    void setReadGroup(String rgId) {
        r1.setAttribute(RG_TAG, rgId);
        r2.setAttribute(RG_TAG, rgId);
    }
    
    @CompileStatic
    int minPos() {
        Math.min(r1?.alignmentStart?:Integer.MAX_VALUE, r2?.alignmentStart?:Integer.MAX_VALUE)
    }
    
    @CompileStatic
    String getReadName() {
        r1 != null ? r1.readName : r2?.readName
    }
    
    /**
     * If one of the pair is missing, return
     * @return
     */
    @CompileStatic
    int getMissingReadPosition() {
        if(r1 != null)
            return r1.mateAlignmentStart
            
        if(r2 != null)
            return r2.mateAlignmentStart
            
        return 0i
    }
    
    /**
     * If one of the pair is missing, return
     * @return
     */
    @CompileStatic
    String getMissingReadReferenceName() {
        if(r1 != null)
            return r1.mateReferenceName
            
        if(r2 != null)
            return r2.mateReferenceName
            
        throw new IllegalStateException("No read populated but read reference name requested")
    }
  
    
    @CompileStatic
    SAMRecord getNonMissingRead() {
        r1 ?: r2
    }
    
    void setR2(SAMRecord r2) {
        assert r2 != null
        this.r2 = r2
    }

    /**
     * Compare this read pair to another read pair by genomic position.
     * <p>
     * Within a chromosome, the R1 position is the primary determinant 
     * for comparision. If R1 positions are equal, R2 position is used,
     * and if these are equal, the read name is compared.
     */
    @CompileStatic
    @Override
    public int compareTo(Object o) {
        
        SAMRecordPair other = (SAMRecordPair)o
        
        int cmp = r1ReferenceIndex - other.r1ReferenceIndex
        if(cmp != 0)
            return cmp
        
        if(r1 != null)
            cmp = this.r1.alignmentStart - other.r1.alignmentStart
            
        if(cmp != 0)
            return cmp
            
        if(r2 == null)
            return 0
            
        int otherAlignmentStart = 0
        if(other.r2) {
            otherAlignmentStart = other.r2.alignmentStart
        }
        
        cmp = this.r2.alignmentStart - otherAlignmentStart
        if(cmp != 0)
            return cmp
            
        return this.readName.compareTo(other.readName)
    }
    
    @CompileStatic
    int getR1ReferenceIndex() {
        r1 != null ? r1.referenceIndex : r2.mateReferenceIndex
    }
    
    @CompileStatic
    int getR2ReferenceIndex() {
        r2 != null ? r2.referenceIndex : r1.mateReferenceIndex
    }
  
    @CompileStatic
    int getR1Pos() {
        r1 != null ? r1.alignmentStart : r2.mateAlignmentStart
    }
  
    @CompileStatic
    int getR2Pos() {
        r2 != null ? r2.alignmentStart : r1.mateAlignmentStart
    }
    
    @CompileStatic
    String getLeftPosition() {
        if(r1ReferenceIndex < r2ReferenceIndex)
            return r1Position
        else
        if(r1ReferenceIndex > r2ReferenceIndex)
            return r2Position
        
        if(r1Pos <= R2Pos)
            r1Position
        else
            r2Position
    }
    
    @CompileStatic
    String getRightPosition() {
        if(r1Pos < R2Pos)
            r2Position
        else
            r1Position
    }
    
    @CompileStatic
    String getR1Position() {
        
        if(r1) {
            r1.referenceName + ':' + r1.alignmentStart
        }
        else {
            r2.mateReferenceName + ':' + r2.mateAlignmentStart
        }
    }
    
    @CompileStatic
    String getR2Position() {
        if(r2) {
            r2.referenceName + ':' + r2.alignmentStart
        }
        else {
            r1.mateReferenceName + ':' + r1.mateAlignmentStart
        }
    }
    
    @CompileStatic
    String toString() {
        this.readName + ": " + leftPosition + " - " + rightPosition
    }
}

