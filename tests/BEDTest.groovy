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


class BEDTest {
    
    String testBed = 
    """
    chr1\t100\t150
    chr2\t190\t250
    chr3\t300\t350
    """.stripIndent().trim()

    @Test
    public void testContains() {
        BED b = new BED(new ByteArrayInputStream(testBed.bytes))
        b.load()
        
        assert b.contains("chr1", 105)
        assert !b.contains("chr1", 99)
        assert b.contains("chr1", 100)
        assert !b.contains("chr1", 150)
        assert !b.contains("chr2", 115)
    }
    
    @Test
    public void testOverlaps() {
        BED b = new BED(new ByteArrayInputStream(testBed2.bytes))
        b.load()
        
        int count = 0
        b.eachOverlap("chr1",115) { String chr, int pos, int end ->
            ++count
        }
        assert count == 2
        
        count = 0 
        b.eachOverlap("chr1",125) { String chr, int pos, int end ->
            ++count
        }
        assert count == 3 
        
        count = 0 
        b.eachOverlap("chr1",1500) { String chr, int pos, int end ->
            ++count
        }
        assert count == 0 
        
        count = 0 
        b.eachOverlap("chr2",125) { String chr, int pos, int end ->
            ++count
        }
        assert count == 0 
    }
    
   String testBed2 = 
    """
    chr1\t100\t350
    chr1\t110\t140
    chr1\t120\t360
    """.stripIndent().trim()
    
    @Test
    void testLongPrecedingRange() {
        BED b = new BED(new ByteArrayInputStream(testBed2.bytes))
        b.load()
        assert b.contains("chr1", 200)
    }
    
    @Test
    void testExtra() {
        BED b = new BED()
        b.add("chr1", 100, 350, "foo")
        b.add("chr1", 110, 250, "bar")
        
        assert b.getExtrasAtPosition("chr1", 200).sort() == ["bar", "foo"]
        
    }
    
    @Test
    void testStartingAt() {
        BED b = new BED()
        b.add("chr1", 100, 350)
         .add("chr1", 110, 250)
       
        assert b.startingAt("chr1", 100) == [100..349] 
        
        b.add("chr1", 100, 370)
         
        assert b.startingAt("chr1", 100) == [100..349, 100..369] 
    }

    @Test
    void testEndingAt() {
        BED b = new BED()
        b.add("chr1", 100, 350)
         .add("chr1", 110, 250)
       
        assert b.endingAt("chr1", 349) == [100..349] 
        
        b.add("chr1", 20, 350)
         
        assert b.endingAt("chr1", 349) == [100..349, 20..349] 
    }
    
    @Test
    void testStartEndSamePos() {
        BED b = new BED()
        b.add("chr1", 100, 350)
         .add("chr1", 350, 500)
         .add("chr1", 100, 120)
         .add("chr1", 350, 450)
       
        assert b.endingAt("chr1", 349) == [100..349] 
        assert b.startingAt("chr1", 100) == [100..349,100..119] 
        assert b.startingAt("chr1", 350) == [350..499,350..449] 
        
//        b.add("chr1", 20, 350)
         
//        assert b.endingAt("chr1", 349) == [100..349, 20..349] 
    }
    
    @Test
    void testWeirdFailingExample() {
        
        def bedFile = """
              chr1\t1276\t1537
              chr1\t860\t1010
              chr1\t915\t1102
              chr1\t1095\t1204
              chr1\t1102\t1383
        """.stripIndent().trim()
                      
        BED b = new BED(new ByteArrayInputStream(bedFile.bytes))
        b.load()
        
        assert b.endingAt("chr1", 1101).size() != 0
    }
    
    @Test
    void testIntersect() {
        BED b = new BED(new ByteArrayInputStream(testBed.bytes))
        b.load()
        
        println b.intersect("chr1", 80, 120)

    }
    
    @Test
    void testSubtractFrom() {
        BED b = new BED()
        b.add("chr1",100,150)
        
        assert b.subtractFrom("chr1", 80, 120) == [80..99]

        b.add("chr1",90,95)
        assert b.subtractFrom("chr1", 80, 120) == [80..89,95..100]
    }
}
