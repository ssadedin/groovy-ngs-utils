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

/**
 * Region of a genome
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class Region extends Expando {
    
    Region() {
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
    
    boolean overlaps(Region other) {
        if(other.chr != this.chr)
            return false
        return this.range.disjoint(other.range)
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
    
    @CompileStatic
    GRange(int from, int to, Object extra) {
        super(from,to)
        this.extra = extra
    }
    
    Object extra    
}
