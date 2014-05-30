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
        assert b.subtractFrom("chr1", 80, 120) == [80..89,95..99]
        
        assert b.subtractFrom("chr1", 92, 120) == [95..99]
    }
    
    @Test
    void testGetOverlaps() {
        BED b = new BED(new ByteArrayInputStream(
          """
          chr1\t100\t150
          chr1\t190\t250
          chr1\t300\t350
          """.stripIndent().trim().bytes
        ))
        b.load()
        
        List<Range> o = b.getOverlaps("chr1", 80,120)
        assert o.size() == 1
        assert o[0] == 100..149
        
        o = b.getOverlaps("chr1", 145,195)
        assert o == [100..149, 190..249]
    } 
    
    @Test
    void testEachUniqueRange() {
        BED b = new BED(new ByteArrayInputStream(
          """
          chr1\t100\t150
          chr1\t190\t250
          chr1\t190\t250
          chr1\t300\t350
          """.stripIndent().trim().bytes
        ))
        b.load()
        
        int count = 0
        Closure c = { chr, start, end ->
            ++count
            println "I was invoked, yes I was: $chr:$start-$end"
        }
        println "Calling with closure ${c.hashCode()}"
        
        b.eachRange(unique:true, c)
        
        assert count == 3
    }
    
    @Test
    void testIterator() {
        BED b = new BED(new ByteArrayInputStream(
          """
          chr1\t100\t150
          chr1\t190\t250
          chr1\t190\t250
          chr1\t300\t350
          """.stripIndent().trim().bytes
        ))
        b.load()
        
        int count = 0
        for(r in b) {
            ++count
        }
        assert b.grep { println "Examining range $it"; it.from == 300 }[0].to == 349
        
        // Conversion to collection is not supposed to be necessary - 
        // I'm not sure why : 'grep' is found but not 'min', both are part of the Iterable
        // interface in Groovy ....
        Iterable i = b as Collection
        assert i.min { it.from }.from == 100
    } 
    
    @Test
    void testReduce() {
        BED b = new BED(new ByteArrayInputStream(
          """
          chr1\t100\t120
          chr1\t140\t210
          chr1\t190\t250
          chr1\t300\t350
          """.stripIndent().trim().bytes
        ))
        b.load()
        
        Regions reduced = b.reduce()
        assert reduced.allRanges["chr1"].size() == 3
    }
    
    @Test
    void testUnique() { 
        BED b = new BED(new ByteArrayInputStream(
            """
            chrX\t8503465\t8503828
            chrX\t8503517\t8503763
            chrX\t8503709\t8503982
            chrX\t8503752\t8503855
            chrX\t8503753\t8503935
            chrX\t8503763\t8504106
            chrX\t8503828\t8504201
            """.stripIndent().trim().bytes)).load()
        
        assert b.endingAt("chrX", 8503981).size() == 1
        
        Regions u = b.unique()
        assert u.endingAt("chrX", 8503981).size() == 1
    }
    
    @Test
    void testIntersectOtherBED() {
      BED b1 = new BED(new ByteArrayInputStream(
          """
          chr1\t100\t120
          chr1\t140\t210
          chr1\t190\t250
          chr1\t300\t350
          """.stripIndent().trim().bytes
        )).load()
        
      BED b2 = new BED(new ByteArrayInputStream(
          """
          chr1\t100\t120
          chr1\t140\t210
          chr1\t190\t250
          chr1\t300\t350
          """.stripIndent().trim().bytes
        )).load()
        
      def b3 = b1.intersect(b2).reduce()
      
      b3.eachRange { println(it.toString()) }
      
      assert b3.iterator().size()>0
        
      // Identical regions intersected should just return the same result
      assert b1.intersect(b2).reduce().size() == b1.reduce().size()
    }
}
