/*
 *  Groovy NGS Utils - Some simple utilites for processing Next Generation Sequencing data.
 *
 *  Copyright (C) 2013 Simon Sadedin, ssadedin<at>gmail.com
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

import groovy.transform.CompileStatic;
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.util.Interval

import java.util.regex.Pattern

@CompileStatic
interface IRegion {
	
	String getChr()
	
	IntRange getRange()
}

/**
 * Region comparator that orders regions first lexically on the contig/chr name
 * then numerically on the positional columns.
 */
@CompileStatic
class RegionComparator implements Comparator<IRegion> {

    @Override
    public int compare(IRegion r1, IRegion r2) {
      if(r1.chr != r2.chr)
        return -r2.chr.compareTo(r1.chr)
      
      if(r1.range.from != r2.range.from)
        return r1.range.from.compareTo(r2.range.from)
      
      return r1.range.to.compareTo(r2.range.to)
    }
    
}

/**
 * Region comparator that orders regions first numerically on the contig/chr name
 * (stripping any non-numeric prefix from `chr` and using special handling for 
 * common chromosomes such as <code>chrX</code>, then numerically on the positional columns).
 */
@CompileStatic
class NumericRegionComparator implements Comparator<IRegion> {

    @Override
    public int compare(IRegion r1, IRegion r2) {
        int result = (int)Math.signum((long)XPos.computePos(r1.chr, r1.range.from) - XPos.computePos(r2.chr, r2.range.from))
        return result
    }
    
}


/**
 * Region of a genome
 * <p>
 * Among all the different region-based objects this is the highest level and most heavy weight
 * representation of a region. The Region is a Groovy {@link groovy.util.Expando} object, so 
 * arbitrary properties can be set on it dynamically.
 * <p>
 * Regions can be created from common string representations of regions:
 * <pre>
 * Region r = new Region("chr1:20000-300000")
 * </pre>
 * Alternatively they can be created directly:
 * <pre>
 * Region r = new Region("chr1", 20000, 300000)
 * </pre>
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class Region extends Expando implements IRegion, Serializable {
    
    public static final long serialVersionUID = 0L
    
    final static Region EMPTY_REGION = new Region("empty", 0..0) {
        @CompileStatic
        long size() {
            return 0l
        }
    }
    
    Region() {
    }
    
    Region(IRegion other) {
        this.chr = other.chr
        GRange newRange = new GRange(other.range.from, other.range.to, this)
        this.range = newRange
    }
    
    Region(Map props=null, String region) {
        parseRegion(region)
        this.range.extra = this
        setProps(props)
    }
    
    Region(Map props=null, String chr, Range range) {
        this.chr = chr
        this.range = range
        setProps(props)
    }
    
    /**
     * Create region on the given chromosome spanning from the given
     * <code>from</code> position to the given <code>to</code> position
     * inclusive.
     */
    Region(Map props=null, String chr, int from, int to) {
        this.chr = chr
        this.range = new GRange(from, to, this)
        setProps(props)
    }
    
    Region(SAMRecord read) {
        this(read.referenceName, read.alignmentStart, read.alignmentEnd)
        this.extra = read
        this.read = read
    }
    
    @CompileStatic
    public void setProps(Map props) {
       if(props == null) 
           return
       for(Map.Entry e in props) {
           this.setProperty((String)e.key, e.value)
       }
    }
    
    @CompileStatic
    void parseRegion(String region) {
        int colonIndex = region.indexOf(":")
        List tokens
        if(colonIndex>0) {
            this.chr = region.substring(0, colonIndex)
            int dashIndex = region.indexOf("-")
            try {
                this.range = new GRange(region.substring(colonIndex+1, dashIndex).toInteger(),region.substring(dashIndex+1).toInteger(), null)
            }
            catch(java.lang.NumberFormatException f) {
                this.range = new GRange(region.substring(colonIndex+1, dashIndex).replace(',','').toInteger(),region.substring(dashIndex+1).replace(',','').toInteger(), null)
            }
        }
        else
        if((tokens = region.tokenize('\t')).size() > 1) {
            this.range = 
                new GRange(tokens[1].toInteger(), 
                           tokens[2].toInteger(), null)
           this.chr = tokens[0]
        }
        else {
           this.chr = region
           this.range = new GRange((int)0,(int)0,null)
        }
    }
    
    // Causes some weird conflict with toString()
    // @Delegate
    IntRange range
    
    String chr
    
    /**
     * The start position of the range represented by this region
     */
    @CompileStatic
    Integer getFrom() {
        return range.from
    }
    
    /**
     * The end position of the range represented by this region
     */
    @CompileStatic
    Integer getTo() {
        return range.to
    }
    
    String startString() {
        "$chr:$from"
    }
    
    String endString() {
        "$chr:$to"
    } 
    
    @CompileStatic
    Region widen(int bp) {
        return widen(bp,bp)
    }
    
    /**
     * @param bp    number of base pairs to widen by
     * 
     * @return  a new Region object representing an interval bp wider than this one,
     *          with the lower boundary being limited to zero
     */
    @CompileStatic
    Region widen(int leftBp, int rightBp) {
        if(range instanceof GRange) {
            return new Region(this.chr, 
                new GRange(Math.max(this.from-leftBp,0), this.to+rightBp, ((GRange)range).extra))
        }
        else {
            Region result = new Region(this.chr, 
                new IntRange(Math.max(this.from-leftBp,0), this.to+rightBp))            
        }
    }
    
    
    /**
     * Return true of this region overlaps the other.
     * <p>
     * Supports the groovy syntax:
     * <pre>
     * Region r1 = new Region("chr1:1-10")
     * Region r2 = new Region("chr1:5-10")
     * 
     * assert r1 in r2
     * </pre>
     * 
     * Note however that any overlap will return true, the region does not 
     * have to be contained entirely.
     * 
     * @param other
     * @return true iff this region overlaps the other
     */
    @CompileStatic
    boolean isCase(IRegion other) {
        return overlaps(other)
    }
    
    /**
     * @return true iff the specified region partially or fully overlaps this one
     */
    @CompileStatic
    boolean overlaps(String chr, int from, int to) {
        if(chr != this.chr)
            return false
        boolean result = this.range.containsWithinBounds(to) || 
               this.range.containsWithinBounds(from) ||
               (this.range.to >= from && this.range.to <= to)
        return result
    }
    
    /**
     * @return true if any region within the given {@link Regions} partially 
     *              or fully overlaps this one
     */
    @CompileStatic
    boolean overlaps(Regions regions) {
        regions.overlaps(this)
    }
    
    /**
     * @return true iff the specified region partially or fully overlaps this one
     */ 
    @CompileStatic
    boolean overlaps(IRegion other) {
        if(other.chr != this.chr)
            return false
        GRange.overlaps(this.range,other.range)
    }
    
    /**
     * @return true iff the specified region partially or fully overlaps this one
     */
    @CompileStatic
    static boolean overlaps(final IRegion r1, final IRegion r2) {
        if(r1.chr != r2.chr)
            return false
        return GRange.overlaps(r1.range,r2.range)
    }
     
    
    /**
     * Compute the fraction of mutual overlap of this region with the other.
     * <p>
     * This is defined as, the minimum fraction of either region's overlap
     * with the other.
     * 
     *               |   20   |
     * |----------------------|
     *           60
     *               |------------------------------|
     *                               80
     * mutual overlap = min(20/60, 20/80) = 0.25
     *                               
     * @param other
     * @return
     */
    @CompileStatic
    double mutualOverlap(IRegion other) {
        
        int ixSize = (int)intersect(other).size()
        
        if(ixSize == 0) 
            return 0.0d
        
        // Since we know that they intersect, the union is just the max coord - min coord
        // Jaccard
//        return ixSize / (Math.max(to, other.range.to) - Math.min(from, other.range.from))
        
        return Math.min(ixSize / (double)size(), ixSize / (double)other.range.size())
    }
    
    /**
     * @return  true iff this region fully encompasses the given region
     */
    @CompileStatic
    boolean spans(IRegion r) {
        (r.chr == this.chr) && (r.range.to <= this.to) && (r.range.from >= this.from)
    }
    
    @CompileStatic
    Region intersect(IRegion other) {
        if(this.chr != other.chr)
            return EMPTY_REGION
            
        if(this.is(EMPTY_REGION))
            return EMPTY_REGION
           
        if(other.is(EMPTY_REGION))
            return EMPTY_REGION
            
        int ixFrom= Math.max(this.from, other.range.from)
        int ixTo = Math.min(this.to, other.range.to)
        if(ixTo<ixFrom)
            return Region.EMPTY_REGION
        return new Region(this.chr, ixFrom..ixTo)
    }
    
    @CompileStatic
    Region union(IRegion other) {
        union(this, other)
    }
 
    Region copy() {
        new Region(chr, range.from..range.to)
    }
    
    @CompileStatic
    Object getExtra() {
        range instanceof GRange ? ((GRange)range).extra : null
    }
    
    /**
     * @return  the number of positions spanned by this region
     */
    @CompileStatic
    long size() {
        range.size()
    }
    
    private static Pattern ALTERNATE_HAPLOTYPE_PATTERN = ~'.*_hap[0-9]*$'
    
    @CompileStatic
    boolean isMinorContig() {
        Region.isMinorContig(this.chr)
    }
    
    /**
     * A hueristic that returns true if this the given contig 
     * represents a non-primary assembly contig or alternate haplotype 
     * in commonly used human genome assemblies.
     */
    @CompileStatic
    static boolean isMinorContig(String chr) {
        chr.startsWith('NC_') ||
        chr.startsWith('GL') ||
        chr.startsWith('Un_') ||
        chr.startsWith('chrUn_') ||
        chr.startsWith('M') ||
        chr.startsWith('chrM') ||
        chr.startsWith('hs3') ||
        chr.endsWith('_random') || 
        chr.endsWith('_alt') ||
        chr.endsWith('_fix') ||
        chr.endsWith('EBV') ||
        chr.startsWith('HLA-') ||
        chr.matches(ALTERNATE_HAPLOTYPE_PATTERN)
    }
    
    @CompileStatic
    int getMidpoint() {
        return (int)((range.to + range.from)/2)
    }
    
    /**
     * There is a design flaw in that because the range is inclusive of both start
     * and end, we cannot represent an empty region. Some functions 
     * therefore return a special designated reigon that signifies "empty".
     * 
     * @return
     */
    boolean isEmpty() {
        this.is(EMPTY_REGION)
    }
    
    @CompileStatic
    Object asType(Class clazz) {
        if(clazz == Interval) {
            return new Interval(chr, range.from, range.to)
        }
        return null
    }
    
    String igv() {
        return "http://localhost:60151/goto?locus=$chr:$from-$to"
    }
    
    @CompileStatic
    Region stripContigPrefix() {
        new Region(this.chr.replace('chr',''), this.range)
    }
    
    String toString() {
        "$chr:$range.from-$range.to"
    }
    
    @CompileStatic
    static Region union(IRegion r1, IRegion r2) {
        assert overlaps(r1,r2)
        new Region(r1.chr, Math.min(r1.range.from, r2.range.from), Math.max(r1.range.to, r2.range.to))
    }
}

/**
 * Range over an interval
 * 
 * @author simon.sadedin@mcri.edu.au
 */
@CompileStatic
class GRange extends IntRange implements Serializable {
    
    public static final long serialVersionUID = 0L
    
    GRange() { // for serialization
        super(0,0)
    }
    
    GRange(int from, int to, Object extra) {
        super(from,to)
        this.extra = extra
    }
    
    @CompileStatic
    boolean spans(IntRange r) {
        return (r.to <= this.to) && (r.from >= this.from)
    }
    
    boolean overlaps(IntRange other) {
        GRange.overlaps(this, other)
    }
    
    @CompileStatic
    GRange intersectRange(IntRange other) {
        GRange r = new GRange(Math.max(this.fromInt, other.fromInt), Math.min(this.toInt, other.toInt), extra)
        if(r.to < r.from)
            return null
        return r
    }
  
    /**
     * @param  other
     * @return a new GRange whose left position is the minimum of the overlap between this
     *         GRange and the other, and whose right position is the maximum of the overlap between
     *         this range and the other.
     */
    @CompileStatic
    GRange intersect(IntRange other) {
        intersectRange(other)
    }
    
    /**
     * @param bp    number of base pairs to widen by
     * 
     * @return  a new GRange object representing an interval bp wider than this one,
     *          with the lower boundary being limited to zero
     */    
    @CompileStatic
    GRange widen(int bp) {
        return widen(bp,bp)
    }
    
    /**
     * @param leftBp    number of base pairs to widen by on the "left" side
     * @param rightBp    number of base pairs to widen by on the "right" side
     * 
     * @return  a new GRange object representing an interval bp wider than this one,
     *          with the lower boundary being limited to zero
     */
    @CompileStatic
    GRange widen(int leftBp, int rightBp) {
        new GRange(Math.max(this.from-leftBp,0), this.to+rightBp, extra)
    }
    
    static boolean overlaps(IntRange a, IntRange b) {
        boolean result = a.containsWithinBounds(b.to) || 
                         a.containsWithinBounds(b.from) ||
                         b.containsWithinBounds(a.to)
        return result
    }
    
    Object extra    
}

class GRegion extends Region {
   GRegion(String chr, Range range) {
        this.chr = chr
        this.range = new GRange(range.from,range.to,this)
   }
}
