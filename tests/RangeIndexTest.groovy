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
    void testOverlapEndPoint() {
        RangeIndex index = new RangeIndex()
        
        [0..50, 20..70, 10..90, 120..130].each {index.add(it)}
         
        assert index.getOverlaps(110,120) == [120..130]
        assert index.getOverlaps(110,119) == []
        
        println index.getOverlaps(130,140) == [120..130]
    }
    
    @Test
    void testOneBaseOverlap() {
        RangeIndex index = new RangeIndex()
        
        [20..20].each {index.add(it)}
        
        index.dump()
        
        assert index.getOverlaps(20,20).size() == 1
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
    void getOverlapsHigherDisjoint() {
       RangeIndex index = new RangeIndex()
       index.ranges[48672868] = [48672868..48672928]
       index.ranges[48672929] = [48672929..48672956]
       
       def overlaps = index.getOverlaps(76763881,76763892)
       assert overlaps.size()==0
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
    
    @Test
    void testRemove() {
       RangeIndex index = new RangeIndex()
        [0..50, 
         70..90, 
         80..85, 
         60..65].each {index.add(it.from, it.to+1)}
            
         index.remove(0..50)
         
         assert index.getOverlaps(30).isEmpty()
         assert index.getOverlaps(82).size() == 2
         index.remove(80..85)
         assert index.getOverlaps(82).size() == 1
         index.remove(70..90)
         assert index.getOverlaps(82).isEmpty()
         
    }
    
    @Test 
    void testIterate() {
       RangeIndex index = new RangeIndex()
       [0..50, 
         70..90, 
         80..85, 
         60..65].each {index.add(it.from, it.to+1)}
            
       for(Range r in index) {
//           println "${r.from} - ${r.to}"
           println r
       }
    }
    
    @Test 
    void testIterateOverlapping() {
       RangeIndex index = new RangeIndex()
       [  
         0..100, 
         50..60, 
         70..80 
       ].each {index.add(it)}
       
       println "Ranges at 0 == ${index.ranges[0]}"
       println "Ranges at 50 == ${index.ranges[50]}"
            
       def ranges = []
       for(Range r in index) {
           ranges.add(r)
       }
       assert ranges.size() == 3
    }
    
    @Test 
    void testReverseIterate() {
       RangeIndex index = new RangeIndex()
       [0..50, 
         70..90, 
         80..85, 
         60..65].each {index.add(it)}
            
       for(Range r in index.reverseIterator()) {
           println "${r.from} - ${r.to}"
       }
    }
    
    @Test 
    void testReverseIterateOverlapping() {
       RangeIndex index = new RangeIndex()
       [  
         0..100, 
         50..60, 
         70..80 
       ].each {index.add(it)}
       
       println "Ranges at 0 == ${index.ranges[0]}"
       println "Ranges at 50 == ${index.ranges[50]}"
            
       def ranges = []
       for(Range r in index.reverseIterator()) {
           println "${r.from} - ${r.to}"
           ranges << r
       }
       
       assert ranges.size() == 3
       assert ranges[0] == 70..80
       assert ranges[1] == 50..60
       assert ranges[2] == 0..100
       
       ranges = index.reverseIteratorAt(50).collect { it }
       assert ranges.size() == 2
       assert ranges[0] == 50..60
       assert ranges[1] == 0..100
       
       ranges = index.reverseIteratorAt(70).collect { it }
       assert ranges.size() == 3
       assert ranges[0] == 70..80
       assert ranges[1] == 50..60
       assert ranges[2] == 0..100
    }
    
    @Test
    void testIteratorAt() {
       RangeIndex index = new RangeIndex()
       [0..50, 
         70..90, 
         80..85, 
         60..65].each {index.add(it)}
            
       def ranges = index.iteratorAt(70).collect { it }
       assert ranges.size == 2
       assert ranges[0] == 70..90
       
       ranges = index.iteratorAt(0).collect { it }
       assert ranges.size == 4
       assert ranges[3] == 80..85
    }
    
    @Test
    void testReverseIteraterAt() {
       RangeIndex index = new RangeIndex()
       [0..50, 
         70..90, 
         80..85, 
         60..65].each {index.add(it)}
            
       def ranges = index.reverseIteratorAt(70).collect { println "$it.from - $it.to"; it }
       assert ranges.size == 3
       assert ranges[0] == 70..90
       assert ranges[1] == 60..65
       assert ranges[2] == 0..50
       
       ranges = index.reverseIteratorAt(0).collect { it }
       assert ranges.size == 1
       assert ranges[0] == 0..50
       
       ranges = index.reverseIteratorAt(60).collect { it }
       assert ranges.size == 2
       assert ranges[0] == 60..65
    }
     
    @Test 
    void testBigIterate() {
//        BED bed = new BED("/Users/simon/work/dsd/batch3/design/tmp.bed")
//        bed.load()
        
       RangeIndex index = new RangeIndex()
        [
          6079..6350,
          6257..6494,
          6923..7296,
          7168..7274,
          7220..7359,
          7226..7401
        ].each {index.add(it.from, it.to+1)}
        
        int count = 0 
        ProgressCounter c = new ProgressCounter()
        for(Range r in index) {
            println "$r.from\t$r.to"
            ++ count
            c.count()
        }
        println "Counted $count ranges"
        assert count == 6
    }
    
    
    @Test 
    void testIterate3() {
//        BED bed = new BED("/Users/simon/work/dsd/batch3/design/tmp.bed")
        BED bed = new BED("/Users/simon/work/dsd/batch3/design/amplicons_with_filtered_enzymes.bed")
        bed.load()
        
        int count = 0 
        ProgressCounter c = new ProgressCounter()
        for(Range r in bed.index["chrX"]) {
//            println "$r.from\t$r.to"
            ++ count
            c.count()
        }
        println "Counted $count ranges"
    }
    
    @Test
    void testRemove2() {
        RangeIndex index = new RangeIndex()
        index.ranges[1581349]=[1581349..1581360]
        index.ranges[1581360]=[1581360..1581360]
        
        index.remove(1581349..1581360)
    }
    
    @Test
    void testRemove3() {
        RangeIndex index = new RangeIndex()
        index.ranges[50155286]=[]
        index.ranges[50155287]=[50155285..50155296]
        index.ranges[50155288]=[50155285..50155296]
        index.ranges[50155289]=[50155285..50155296]
        index.ranges[50155290]=[50155285..50155296]
        index.ranges[50155291]=[50155285..50155296]
        index.ranges[50155292]=[50155285..50155296]
        index.ranges[50155293]=[50155285..50155296]
        index.ranges[50155294]=[50155285..50155296]
        index.ranges[50155295]=[50155285..50155296]
        
        index.remove(50155285..50155296)
    }
    
    @Test
    void getOverlaps4() {
        RangeIndex index = new RangeIndex()
        index.ranges[1273475] = [1273475..1273503]
        index.ranges[1273504]= []
        
        assert index.getOverlaps(1273503,1273504).size() == 1
    }
    
    @Test
    void testIntersection() {
        RangeIndex index = new RangeIndex()
        [0..50, 20..70, 10..90, 120..130].each {index.add(it)}
        
        assert index.intersect(100,125) == [120..125]
        assert index.intersect(130,135) == [130..130]
        assert index.intersect(80,125) == [80..90,120..125]
    }
    
    @Test
    void testGetOverlaps() {
        RangeIndex index = new RangeIndex()
            [ 
             32503143..32503704,
             45981437..45981700
            ].each { index.add(it) }
            
        assert index.getOverlaps(32503143, 32503704).size() > 0
    }
    
    // 1273475:[1273475..1273503], 1273503:[1273503..1273503], 1273504:[], 1284266:[1284266..1284268], 1284268:[1284268..1284268], 1284269:[], 1575676:[1575676..1575676], 1575677:[], 1581349:[1581349..1581360], 1581360:[1581360..1581360], 1581361:[], 182992872:[182992872..182992971], 182992971:[182992971..182992971], 182992972:[]]:

    
}
