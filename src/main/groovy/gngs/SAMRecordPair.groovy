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
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SAMTagUtil
import htsjdk.samtools.SAMUtils
import htsjdk.samtools.util.SequenceUtil

@CompileStatic
interface ReadPair {
    boolean isChimeric()
    boolean getUnmapped()
    boolean notInRegions(Regions regions)
}

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
@CompileStatic
class SAMRecordPair implements Comparable, ReadPair {
    
    public SAMRecord r1
    
    public SAMRecord r2
    
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
    
    @CompileStatic
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
    String getR1ReferenceName() {
        r1 != null ? r1.referenceName : r2.mateReferenceName
    } 
    
    @CompileStatic
    String getR2ReferenceName() {
        r2 != null ? r2.referenceName : r1.mateReferenceName
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
        
        if(r1 != null) {
            r1.referenceName + ':' + r1.alignmentStart
        }
        else {
            r2.mateReferenceName + ':' + r2.mateAlignmentStart
        }
    }
    
    @CompileStatic
    String getR2Position() {
        if(r2 != null) {
            r2.referenceName + ':' + r2.alignmentStart
        }
        else {
            r1.mateReferenceName + ':' + r1.mateAlignmentStart
        }
    }
    
    @CompileStatic
    String toString() {
        if(this.r1ReferenceIndex == this.r2ReferenceIndex)
            this.readName + " @ " + r1Position + "-" + r2Pos
        else
            this.readName + " @ " + leftPosition + "-" + rightPosition
    }
    
    @CompileStatic
    int getReadLength() {
        if(r1 != null)
            r1.readLength
        else
            r2.readLength
    }
    
    /**
     * @return true if either r1 or r2 is unmapped
     */
    @CompileStatic
    boolean getUnmapped() {
        if(r1 != null)
            return r1.readUnmappedFlag || r1.mateUnmappedFlag
        else
            return r2.readUnmappedFlag || r2.mateUnmappedFlag
    }
    
    @CompileStatic
    boolean overlaps(Region region) {
        String chr = this.r1ReferenceName
        int r1p = r1Pos
        int len = this.readLength
        if(region.overlaps(chr, r1p, r1p+len)) 
            return true
        
        int r2p = r2Pos
        return region.overlaps(r2ReferenceName, r2p, r2p+len)
    }
    
    boolean spanOverlaps(Region region) {
        int r1p = r1Pos
        int r2p = r2Pos
        
        String chr = this.r1ReferenceName
        if(this.isChimeric()) {
            if(r1p < r2p)
                return (chr == region.chr) && (r1p > region.from) && (r1p < region.to)
            else
                return (chr == region.chr) && (r2p > region.from) && (r1p < region.to)
        }
        else {
            if(r1p < r2p)
                return region.overlaps(chr, r1p, r2p)
            else
                return region.overlaps(chr, r2p, r1p)
        }
    }
    
    @CompileStatic
    void write(Writer out) {
        
        StringBuilder b = new StringBuilder(500)
        b.setLength(0)
        this.appendTo(b,b)
        out.println(b.toString())
    }
    
    @CompileStatic
    void appendTo(final StringBuilder b1, StringBuilder b2, boolean addPosition=false) {
        String r1Name = r1.readName
        StringBuilder nameSuffix = null
        
        if(addPosition) { 
            nameSuffix = new StringBuilder()
            nameSuffix.append(':')
            nameSuffix.append(r1.referenceName)
            nameSuffix.append(':')
            nameSuffix.append(r1.alignmentStart)
            nameSuffix.append(':')
            if(r1.referenceIndex != r2.referenceIndex)
                nameSuffix.append(r2.referenceName)
            nameSuffix.append(r2.alignmentStart)
        }
           
        // R1
        b1.append()
        b1.append('@')
        b1.append(r1Name)
        if(nameSuffix != null) { 
            b1.append((CharSequence)nameSuffix)
        }        
        b1.append(' 1:N:0:1\n')
        b1.append(r1.readString)
        b1.append('\n+\n')
        b1.append(r1.baseQualityString)
        b1.append('\n')
            
        // R2
        b2.append('@')
        b2.append(r1Name)
        if(nameSuffix != null) { 
            b2.append((CharSequence)nameSuffix)
        }        
        
        b2.append(' 2:N:0:1\n')
        b2.append(SequenceUtil.reverseComplement(r2.readString))
        b2.append('\n+\n')
        
        byte [] bq = r2.baseQualities
        for(int i=bq.length-1; i>=0; --i) {
            b2.append(SAMUtils.phredToFastq(bq[i]))
        }
        b2.append('\n')
    }

    @CompileStatic
    @Override
    public boolean notInRegions(Regions regions) {
        
        final String chr1 = this.r1ReferenceName
        final int r1p = this.r1.alignmentStart
        final int len = this.readLength
        if(regions.overlaps(chr1, r1p , r1p+len))
            return false
        
        final String chr2 = this.r2ReferenceName?:chr1
        final int r2p = this.r2Pos
        if(regions.overlaps(chr2, r2p , r2p+len)) {
            return false 
        }
        
        return true 
    }
}

