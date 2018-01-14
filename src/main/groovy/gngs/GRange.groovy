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
import java.util.regex.Pattern

interface IRegion {
	
	String getChr()
	
	IntRange getRange()
}

@CompileStatic
class RegionComparator implements Comparator<Region> {

    @Override
    public int compare(Region r1, Region r2) {
      if(r1.chr != r2.chr)
        return -r2.chr.compareTo(r1.chr)
      
      if(r1.from != r2.from)
        return r1.from.compareTo(r2.from)
      
      return r1.to.compareTo(r2.to)
    }
    
}

@CompileStatic
class NumericRegionComparator implements Comparator<Region> {

    @Override
    public int compare(Region r1, Region r2) {
        int result = (int)Math.signum((long)XPos.computePos(r1.chr, r1.from) - XPos.computePos(r2.chr, r2.from))
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
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class Region extends Expando implements IRegion, Serializable {
    
    public static final long serialVersionUID = 0L
    
    final static Region EMPTY_REGION = new Region("empty", 0..0)
    
    Region() {
    }
    
    Region(String region) {
        parseRegion(region)
        this.range.extra = this
    }
    
    Region(String chr, Range range) {
        this.chr = chr
        this.range = range
    }
    
    Region(String chr, int from, int to) {
        this.chr = chr
        this.range = (from..to)
    }
  
    
    @CompileStatic
    void parseRegion(String region) {
        int colonIndex = region.indexOf(":")
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
        else {
           this.chr = region
           this.range = new GRange((int)0,(int)0,null)
        }
    }
    
    // Causes some weird conflict with toString()
    // @Delegate
    IntRange range
    
    String chr
    
    @CompileStatic
    Integer getFrom() {
        return range.from
    }
    
    @CompileStatic
    Integer getTo() {
        return range.to
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
     * @return  true iff this region fully encompasses the given region
     */
    @CompileStatic
    boolean spans(IRegion r) {
        (r.chr == this.chr) && (r.range.to <= this.to) && (r.range.from >= this.from)
    }
    
    Region intersect(IRegion other) {
        if(this.chr != other.chr)
            return EMPTY_REGION
            
        Region r = new Region(this.chr, Math.max(this.from, other.from)..Math.min(this.to, other.to))
        if(r.to < r.from)
            return Region.EMPTY_REGION
        return r
    }
 
    Region copy() {
        new Region(chr, range.from..range.to)
    }
    
    Object getExtra() {
        range instanceof GRange ? range.extra : null
    }
    
    /**
     * @return  the number of positions spanned by this region
     */
    @CompileStatic
    long size() {
        range.size()
    }
    
    private static Pattern ALTERNATE_HAPLOTYPE_PATTERN = ~'.*_hap[0-9]*$'
    
    boolean isMinorContig() {
        Region.isMinorContig(this.chr)
    }
    
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
        chr.matches(ALTERNATE_HAPLOTYPE_PATTERN)
    }
    
    String toString() {
        "$chr:$range.from-$range.to"
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
    
    boolean spans(IntRange r) {
        return (r.to <= this.to) && (r.from >= this.from)
    }
    
    boolean overlaps(IntRange other) {
        GRange.overlaps(this, other)
    }
    
    GRange intersect(IntRange other) {
        GRange r = new GRange(Math.max(this.fromInt, other.fromInt), Math.min(this.toInt, other.toInt), extra)
        if(r.to < r.from)
            return null
        return r
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
