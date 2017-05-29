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
        return r1.chr.compareTo(r2.chr)
      
      if(r1.from != r2.from)
        return r1.from.compareTo(r2.from)
      
      return r1.to.compareTo(r2.to)
    }
    
}


/**
 * Region of a genome
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
    
    Integer getFrom() {
        return range.from
    }
    
    Integer getTo() {
        return range.to
    }
    
    boolean isCase(IRegion other) {
        return overlaps(other)
    }
    
    @CompileStatic
    boolean overlaps(String chr, int from, int to) {
        if(chr != this.chr)
            return false
        boolean result = this.range.containsWithinBounds(to) || 
               this.range.containsWithinBounds(from) ||
               (this.range.to >= from && this.range.to <= to)
        return result
    }
    
    @CompileStatic
    boolean overlaps(Regions regions) {
        regions.overlaps(this)
    }
    
    @CompileStatic
    boolean overlaps(IRegion other) {
        if(other.chr != this.chr)
            return false
        GRange.overlaps(this.range,other.range)
    }
    
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
class GRange extends IntRange implements Serializable {
    
    public static final long serialVersionUID = 0L
    
    GRange() { // for serialization
        super("chrN",0,0)
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
        GRange r = new GRange(Math.max(this.from, other.from), Math.min(this.to, other.to), extra)
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
