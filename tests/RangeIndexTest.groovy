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
    
    @Test
    void testNextRange() {
        RangeIndex index = new RangeIndex()
        [0..50, 
         70..90, 
         60..65].each {index.add(it.from, it.to)}
         
         println index.dump()
        
        assert index.nextRange(60).from == 70
        assert index.nextRange(25).from == 60
        assert index.nextRange(69).from == 70
        
    }
    
    @Test
    void testPreviousRange() {
        RangeIndex index = new RangeIndex()
        [0..50, 
         70..90, 
         60..65].each {index.add(it.from, it.to)}
         
         println index.dump()
        
        assert index.previousRange(59).from == 0
        assert index.previousRange(70).from == 60
        assert index.previousRange(69).from == 60
        assert index.previousRange(65).from == 60
        assert index.previousRange(64).from == 60
        assert index.previousRange(63).from == 0
        
    } 
    
    @Test
    void testNearest() {
       RangeIndex index = new RangeIndex()
        [0..50, 
         70..90, 
         80..85, 
         60..65].each {index.add(it.from, it.to)}
            
        assert index.nearest(51).from == 0
        assert index.nearest(30).from == 0
        assert index.nearest(54).from == 0
        assert index.nearest(82).from == 80
    }
}
