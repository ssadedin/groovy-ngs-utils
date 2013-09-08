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

import static org.junit.Assert.*;

import org.junit.Test;

class RangeIndexTest {
    
    @Test
    void testSimpleOverlap() {
        
        RangeIndex index = new RangeIndex()
        
        [0..50, 20..70, 10..90].each {index.add(it)}
        
        index.dump()
        
        assert index.getOverlaps(5).size() == 1
        assert index.getOverlaps(15).size() == 2
        assert index.getOverlaps(25).size() == 3
        
        assert 5 in index
        assert !(360 in index)
        
        // 3 breakpoints, so should be 3 entries
//        assert index.ranges.size() == 3
    }
    
    @Test
    void testDisjointRange() {
        RangeIndex index = new RangeIndex()
        
        [0..50, 70..90].each {index.add(it)}
        
        index.dump()
        
        assert 20 in index
        assert !(60 in index)
        
        // 3 breakpoints, so should be 3 entries
//        assert index.ranges.size() == 3        
    }
    
    @Test
    void testBoundaries() {
        RangeIndex index = new RangeIndex()
        
        [0..50, 70..90].each {index.add(it.from, it.to)}
        
        index.dump()
        
        assert 0 in index
        assert !(50 in index)
    }
    
    @Test
    void testEqualStart() {
        RangeIndex index = new RangeIndex()
        
        [0..50, 0..90].each {index.add(it.from, it.to)}
        
        index.dump()
        
        assert 0 in index
        assert 50 in index
        assert 70 in index
        assert !(91 in index)
        assert !(90 in index)
        
        assert index.getOverlaps(20).size() == 2
        assert index.getOverlaps(70).size() == 1
        assert index.getOverlaps(90).size() == 0
    } 
    
    @Test
    void testFillGap() {
        RangeIndex index = new RangeIndex()
        [0..50, 70..90, 60..65].each {index.add(it.from, it.to)}
    }
}
