package gngs
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

import gngs.BED
import gngs.ProgressCounter
import gngs.RangeIndex
import gngs.IntervalTreeRangeIndex
import gngs.IRangeIndex

class IRangeIndexTest {
    
    /**
     * Helper method to create the IRangeIndex implementation we want to test
     */
    IRangeIndex createIndex() {
        // Change this line to test different implementations
        return new IntervalTreeRangeIndex()
    }
    
    @Test
    void testSimpleOverlap() {
        println "Running test ..."
        
        IRangeIndex index = createIndex()
        
        [0..50, 20..70, 10..90].each {index.add(it)}
        
        assert index.getOverlaps(5).size() == 1
        assert index.getOverlaps(15).size() == 2
        assert index.getOverlaps(25).size() == 3
        
        assert 5 in index
        assert !(360 in index)
    }
    
    @Test
    void testOverlapEndPoint() {
        IRangeIndex index = createIndex()
        
        [0..50, 20..70, 10..90, 120..130].each {index.add(it)}
         
        assert index.getOverlaps(110,120) == [120..130]
        assert index.getOverlaps(110,119) == []
        
        println index.getOverlaps(130,140) == [120..130]
    }
    
    @Test
    void testOneBaseOverlap() {
        IRangeIndex index = createIndex()
        
        [20..20].each {index.add(it)}
        
        assert index.getOverlaps(20,20).size() == 1
    }
    
    @Test
    void testDisjointRange() {
        IRangeIndex index = createIndex()
        
        [0..50, 70..90].each {index.add(it)}
        
        assert 20 in index
        assert !(60 in index)
    }
    
    @Test
    void testBoundaries() {
        IRangeIndex index = createIndex()
        
        [0..50, 70..90].each {index.add(it.from, it.to)}
        
        assert 0 in index
        assert !(50 in index)
    }
    
    @Test
    void testEqualStart() {
        IRangeIndex index = createIndex()
        
        [0..50, 0..90].each {index.add(it.from, it.to)}
        
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
        IRangeIndex index = createIndex()
        [0..50, 70..90, 60..65].each {index.add(it.from, it.to)}
    }
    
    @Test
    void testNextRange() {
        IRangeIndex index = createIndex()
        [0..50, 
         70..90, 
         60..65].each {index.add(it.from, it.to)}
        
        assert index.nextRange(60).from == 70
        assert index.nextRange(25).from == 60
        assert index.nextRange(69).from == 70
    }
    
    @Test
    void testPreviousRange() {
        IRangeIndex index = createIndex()
        [0..50, 
         70..90, 
         60..65].each {index.add(it.from, it.to)}
        
        assert index.previousRange(59).from == 0
        assert index.previousRange(70).from == 60
        assert index.previousRange(69).from == 60
        assert index.previousRange(65).from == 60
        assert index.previousRange(64).from == 60
        assert index.previousRange(63).from == 0
    } 
    
    @Test
    void getOverlapsHigherDisjoint() {
       IRangeIndex index = createIndex()
       index.add(48672868..48672928)
       index.add(48672929..48672956)
       
       def overlaps = index.getOverlaps(76763881,76763892)
       assert overlaps.size()==0
    }
    
//    @Test
    void testNearest() {
       IRangeIndex index = createIndex()
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
       IRangeIndex index = createIndex()
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
       IRangeIndex index = createIndex()
       [0..50, 
         70..90, 
         80..85, 
         60..65].each {index.add(it.from, it.to+1)}
            
       for(Range r in index) {
           println r
       }
    }
    
    @Test 
    void testIterateOverlapping() {
       IRangeIndex index = createIndex()
       [  
         0..100, 
         50..60, 
         70..80 
       ].each {index.add(it)}
            
       def ranges = []
       for(Range r in index) {
           ranges.add(r)
       }
       assert ranges.size() == 3
    }
    
//    @Test 
    void testReverseIterate() {
       IRangeIndex index = createIndex()
       [0..50, 
         70..90, 
         80..85, 
         60..65].each {index.add(it)}
            
       for(Range r in index.reverseIterator()) {
           println "${r.from} - ${r.to}"
       }
    }
    
//    @Test 
    void testReverseIterateOverlapping() {
       IRangeIndex index = createIndex()
       [  
         0..100, 
         50..60, 
         70..80 
       ].each {index.add(it)}
            
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
    
//    @Test
    void testIteratorAt() {
       IRangeIndex index = createIndex()
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
    
//    @Test
    void testReverseIteraterAt() {
       IRangeIndex index = createIndex()
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
       IRangeIndex index = createIndex()
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
    void testIntersection() {
        IRangeIndex index = createIndex()
        [0..50, 20..70, 10..90, 120..130].each {index.add(it)}
        
        assert index.intersect(100,125) == [120..125]
        assert index.intersect(130,135) == [130..130]
        assert index.intersect(80,125) == [80..90,120..125]
    }
    
    @Test
    void testGetOverlaps() {
        IRangeIndex index = createIndex()
            [ 
             32503143..32503704,
             45981437..45981700
            ].each { index.add(it) }
            
        assert index.getOverlaps(32503143, 32503704).size() > 0
    }
    
    @Test
    void testSubtractFrom() {
        IRangeIndex index = createIndex()
        [
            0..100
        ].each { index.add(it) }
         
        List<Range> result = index.subtractFrom(0, 100)
       
        assert result.size() == 0
        
        result = index.subtractFrom(0, 200)
        assert result.size() == 1
        assert result[0].from == 101
        
        result = index.subtractFrom(0, 200)
        assert result.size() == 1
        assert result[0].from == 101
        
        index.add(150..160)

        result = index.subtractFrom(0, 200)
        assert result.size() == 2
        assert result[0].from == 101
        assert result[0].to == 149
        assert result[1].from == 161
        assert result[1].to == 200
    }
    
    @Test 
    void testReduce() {
        IRangeIndex index = createIndex()
        
        // Add overlapping ranges
        index.add(0..10)
        index.add(5..15) 
        index.add(20..30)
        index.add(25..35)
        index.add(50..60)
        
        // Reduce the ranges
        IRangeIndex reduced = index.reduce()
        
        // Should combine to 3 ranges:
        // 0-15 (from merging 0-10 and 5-15)
        // 20-35 (from merging 20-30 and 25-35) 
        // 50-60 (unchanged)
        List<IntRange> ranges = reduced.collect { it }
        
        assert ranges.size() == 3
        assert ranges[0] == 0..15
        assert ranges[1] == 20..35  
        assert ranges[2] == 50..60
        
        // Test with GRanges and reducer
        index = createIndex()
        index.add(0, 10, "A")
        index.add(5, 15, "B") // note: this method subtracts 1 from end point, consistent with original RangeIndex
        
        reduced = index.reduce { r1, r2 -> r1.extra + r2.extra }
        ranges = reduced.collect { it }
        
        assert ranges.size() == 1
        assert ranges[0].from == 0
        assert ranges[0].to == 14
        assert ranges[0] instanceof GRange
        assert ((GRange)ranges[0]).extra == "AB"
    }
    
//    @Test
    void testCoverage() {
       IRangeIndex index = createIndex()
        [0..50, 
         70..90, 
         80..85, 
         60..65].each {index.add(it)}
         
        List cov = index.coverage()
        
        cov.each { println it.from + " - " + it.to + " : " + it.extra }
    }
    
    @Test
    void testOverlapsBug() {
        
        BED clipped = new BED("tests/data/small.overlaps.bed").load()
        
        println "Source regions:"
        clipped.each { println it.from + '-' + it.to }
        
        List overlaps = clipped.index["chr22"].getOverlaps(26204614,26205124)
        
        println "Overlaps are: " 
        overlaps.each { println it.from + '-' + it.to }
        
        assert overlaps.size() == 3
    }
    
    @Test
    void testGetOverlapsDupes() {
       IRangeIndex index = createIndex()
        [0..50, 
         70..90, 
         70..90, 
         160..165].each {index.add(it)}
         
       List overlaps = index.getOverlaps(60,100)
       
       assert overlaps.size() == 2
    }
    
    @Test
    void testAdjacentRegionOverlaps() {
        
        IRangeIndex index = createIndex()
        [
            2..3,
            10..20,
            5..10
        ]
        .each {index.add(it)}
   
        assert index.getOverlaps(12,15).size()>0
    }
    
//    @Test
    void testDistanceTo() {
        IRangeIndex index = createIndex()
        [
            10..20,
            30..40,
            50..60
        ]
        .each {index.add(it)}
         
        assert index.distanceTo(55) == 0
        assert index.distanceTo(50) == 0
        assert index.distanceTo(63) == 3
        assert index.distanceTo(48) == 2
        assert index.distanceTo(5) == 5
        assert index.distanceTo(0) == 10
    }
}
