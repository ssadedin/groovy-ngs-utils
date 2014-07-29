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

interface IRegion {
	
	String getChr()
	
	IntRange getRange()
}


/**
 * Region of a genome
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class Region extends Expando implements IRegion {
    
    Region() {
    }
    
    Region(String region) {
        int colonIndex = region.indexOf(":")
        this.chr = region.substring(0, colonIndex)
        int dashIndex = region.indexOf("-")
        this.range = region.substring(colonIndex+1, dashIndex).toInteger()..region.substring(dashIndex+1).toInteger()
    }
    
    Region(String chr, Range range) {
        this.chr = chr
        this.range = range
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
    boolean overlaps(IRegion other) {
        if(other.chr != this.chr)
            return false
        boolean result = this.range.containsWithinBounds(other.range.to) || 
               this.range.containsWithinBounds(other.range.from) ||
               other.range.containsWithinBounds(this.to)
        return result
    }
    
    boolean spans(IRegion r) {
        (r.chr == this.chr) && (r.range.to <= this.to) && (r.range.from >= this.from)
    }
    
    Region copy() {
        new Region(chr, range.from..range.to)
    }
    
    Object getExtra() {
        range instanceof GRange ? range.extra : null
    }
    
    int size() {
        range.size()
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
class GRange extends IntRange {
    
    GRange(int from, int to, Object extra) {
        super(from,to)
        this.extra = extra
    }
    
    boolean spans(IntRange r) {
        return (r.to <= this.to) && (r.from >= this.from)
    }
    
    Object extra    
}
