import static org.junit.Assert.*;

import org.junit.Test;

import gngs.*

/**
 * Note: most of the tests for Regions are in fact in the {@link BEDTest}
 * test suite because all the original functionality originally was
 * developed as part of the {@link BED} class.
 * 
 * @author simon.sadedin@mcri.edu.au
 *
 */
class RegionsTest {

    @Test
    void testSubtract() {
        
        Regions regions = new Regions([
            new Region("chr1",100..140),   
            new Region("chr1",120..160)
        ])
        
        Regions other = new Regions([
            new Region("chr1",130..150)
        ])
         
        Regions result = regions.subtract(other)
        
        result.each { println(it) }
        
        result.each { println "from: " + it.from + "-" + it.to }
        
        assert result.find { it.from == 100 && it.to == 129 }
        assert result.find { it.from == 120 && it.to == 129 }
        assert result.find { it.from == 151 && it.to == 160 }
        
    }
    
    @Test
    void testSubtractSame() {
        Regions regions = new Regions([
            new Region("chr1",1..5)
        ])
        
        Regions other = new Regions([
            new Region("chr1",1..5)
        ])
         
        Regions result = regions.subtract(other)
        assert result.size() == 0
    }
    
    @Test
    void testSubtract2Same() {
        Regions regions = new Regions([
            new Region("chr1",100..200),
            new Region("chr1",300..500)
        ])
        
        Regions other = new Regions([
            new Region("chr1",100..200),
            new Region("chr1",300..500)
        ])
         
        Regions result = regions.subtract(other)
        assert result.size() == 0
    }
    
    @Test
    void testSubtract1bpDiff() {
        Regions regions = new Regions([
            new Region("chr1",1..6)
        ])
        
        Regions other = new Regions([
            new Region("chr1",1..5)
        ])
         
        Regions result = regions.subtract(other)
        assert result.size() == 1
        assert result[0].size() == 1
    }  
    
    @Test
    void 'test that subtracting 1bp regions removes the regions'() {
        Regions regions = new Regions([
            new Region("chr1",5..5),
            new Region("chr1",8..8)
        ])
        
        Regions other = new Regions([
            new Region("chr1",5..5),
            new Region("chr1",8..8)
        ])
         
        Regions result = regions.subtract(other)
        assert result.size() == 0
        assert result.numberOfRanges == 0
    }  
    
    
//    @Test
    void testInfiniteLoop() {
        
        Regions amplicons = new BED(withExtra:true, "/Users/simon/work/dsd/batch3/simulation/design/amplicons_with_enzymes.bed").load()
        Regions fp = new BED(withExtra:true, "/Users/simon/work/dsd/batch3/simulation/sim3/sim3.angel.cnvs.roc.bed.0.fp.bed").load()
        def inds = amplicons.findIndexValues { fp[1].overlaps(it) }
        
        println inds.size()
    }
    
    @Test
    void testReduceOverlap() {
        
        String testRegions = """
            chr1:32503143-32503704
            chr1:45981437-45981700
            chr1:154918251-154918854
            chr10:161276784-161277177
            chr11:103435986-103436379
            chr11:534218-534432
            chr11:63767589-63767885
            chr11:65375513-65375697
            chr11:65427598-65427826
            chr12:66391794-66392239
        """.trim().stripIndent()
                
        Regions regions = new Regions()
        testRegions.trim().eachLine {
            regions.addRegion(new Region(it.trim()))
        }
        
//        println "All regions are: " + regions*.toString()
        
        Regions reduced = regions.reduce()
        
//        println "Reduced regions are " + reduced*.toString()
        for(Region r in reduced) {
            println "Overlaps of $r are " + regions.getOverlaps(r)
            assert regions.getOverlaps(r).size() > 0
        }
    }
    
    @Test
    void testOverlapBug() {
        
        def testRegions = [
            "chrX:1426-2499",  // first region
            "chrX:2225-3883",  // overlaps region 1
            "chrX:3550-4587",  // overlaps region 2
            ]
        
        Regions regions = new Regions()
        testRegions.each { regions.addRegion(new Region(it)) }
        
        Regions reduced = regions.reduce()
        
        println reduced*.toString()
        
        assert reduced.numberOfRanges == 1
        assert reduced[0].to == 4587
        
        Region target = new Region("chrX:3272-4352")
        println reduced.getOverlaps(target)
    }
    
    @Test
    void testAddRegionProperties() {
        Regions regions = new Regions();
        Region r = new Region("chr1",1..100)
        r.foo = "bar"
        regions.addRegion(r)
        
        // We want the property we set on the region to be available when we query by region
        assert regions[0].foo == "bar"
    }
    
    @Test
    void testPreserveRegionProperties() {
        Regions regions = new Regions();
        Region region = new Region("chr1",1..100)
        regions.addRegion(region)
        
        regions.each { Region r ->
            r.foo = "bar"
        }
        
        regions.each { Region r ->
            assert r.foo == "bar"
        }
    }
    
    @Test
    void testLastSmallRegion() {
        Regions r = new Regions();
        [
            "chr1:198608098-198726606",
            "chr1:198608098-198664301"
        ].each { r.addRegion(new Region(it)) }
        
        Region ptprc = new Region("chr1:198700435-198725412")
        
        assert !r.getOverlaps(ptprc).empty
    }
    
    @Test
    void testWindow() {
        
        Regions regions = new Regions([
            new Region("chr1",100..140),   
            new Region("chr1",120..160),
            new Region("chr1",180..190),
            new Region("chr1",200..220),
            new Region("chr1",240..260),
            new Region("chr1",280..290)
        ])
        
        def window = regions.window(regions[2], 2) as Regions
        
        println window.join(',')
        
        assert window[0].from == 100
        assert window.numberOfRanges == 5
        assert window[-1].from == 240
    }
    
    @Test
    void testSort() {
        Regions r = 
            new Regions([new Region("chr1:100-200"), 
                         new Region("chr1:50-70"), 
                         new Region("chr1:100-300"), 
                         new Region("chr2:100-400"),
                         new Region("chr1:500-700") 
                         ])
        
       Regions sorted = r.toSorted(new RegionComparator())
       
       assert sorted[0].from == 50
       assert sorted[0].chr == "chr1"
       assert sorted[-1].chr == "chr2"
    } 
    
    @Test
    void testCoverage() {
        Regions r = 
            new Regions([new Region("chr1:100-200"), 
                         new Region("chr1:50-70"), 
                         new Region("chr1:100-300"), 
                         new Region("chr2:100-400"),
                         new Region("chr1:500-700") 
                         ])
         
       Regions cov = r.coverage()
       
       assert cov.find { it.overlaps('chr1', 60, 61) }.extra == 1 
       assert cov.find { it.overlaps('chr1', 110, 115) }.extra == 2 
    }
    
// This behavior unfortunately forces the Regions to not set 
// the added region as the extra info for the Region object,
// which in turn means that iteration, etc, do not return the original
// regions, which is very confusing. Therefore the behavior is changing such
// that all regions added to a Regions object have themselves set as the extra
// object on their ranges.
//    
//    @Test
//    void 'extra should not be lost by reduce'() {
//        
//      Region r1 = new Region('chr1',new GRange(1,10,'foo'))
//      assert r1.extra == 'foo'
//      
//      Regions regions = new Regions([
//        r1,
//        new Region('chr1', new GRange(5,25,'bar'))
//      ])
//                             
//      Regions red = regions.reduce()
//      red.each { r ->  println "$r \t: $r.extra" }
//      assert red[0].extra == 'foo'
//      
//      Regions red2 = regions.reduce { e1, e2 ->
//          e2.extra
//      }
//      
//      assert red2[0].extra == 'bar'
//      
//    }
    
    
    @Test
    void 'properties should not be lost by reduce'() {
        
      Region r1 = new Region('chr1',1..10,foo:1)
      assert r1.foo == 1
      
      Regions regions = new Regions([
        r1,
        new Region('chr1', 5..25, bar:2)
      ])
                             
      Regions red = regions.reduce()
      red.each { r ->  println "$r \t: $r.foo" }
      
      assert red[0].foo == 1
      
      Regions red2 = regions.reduce { e1, e2 ->
          e2.extra
      }
      
      assert red2[0].bar == 2
      
    } 
    
    @Test
    void testForward() {
        Regions regions = new Regions([
            new Region("chr1",100..120),   
            new Region("chr1",130..160),
            new Region("chr1",190..210),
            new Region("chr1",215..220),
            new Region("chr1",230..260),
            new Region("chr1",270..290)
        ])
        
        
        // One region
        assert regions.forward("chr1", 212, 1) == 215..220
        assert regions.backward("chr1", 212, 1) == 190..210
        
        // Two regions
        assert regions.forward("chr1", 212, 2) == 230..260
        assert regions.backward("chr1", 212, 2) == 130..160
        
        // From middle of region
        assert regions.forward("chr1", 200, 1) == 215..220
        assert regions.backward("chr1", 200, 1) == 190..210
         
    }
    
//    @Test
    void testWiden() {
        BED bed = new BED("/Users/simon.sadedin/work/stretch/chimerics/chimeric_regions.bed").load()
        
        Utils.time("Widen bed file by 100bp") {
            bed.widen(100)
        }
    }
    
    @Test
    void addExistingRegions() {
        Region foo = new Region("chr1", 1..10, foo: 1)
        
        
//        Regions regions = new Regions()
//        regions.addRegion(foo)
        
        Regions regions = [foo] as Regions
        assert regions.startingAt("chr1", 1)[0].extra == foo
    }
    
    void 'enhanced regions should retain properties during multiple iterations'() {
        BED regions = new Regions()
        regions.add('chr1', 8, 20)
        regions = regions.enhance()
        
        for(Region region in regions) {
            region.foo = 'bar'
        }
        
//        for(Region region in regions) {
//            assert region.foo == 'bar'
//        }
        
        List<Range> overlaps = regions.getOverlaps('chr1',10,12)
        
        assert overlaps.size()==1
        assert overlaps[0].extra != null
        
        Region r = overlaps[0].extra
        
        assert r.foo == 'bar'
        
    }
    
    @Test
    void testWidenPreservesExtra() {
        
        Region r1 = new Region('chr1', 1, 100, foo: "hello")
        Region r2 = new Region('chr1', 50, 120, foo: "world")
        
        Regions rr = [r1,r2] as Regions
        
        Regions w = rr.widen(10)
        
        assert w[0].extra.foo == "hello"
        
        Regions wr = rr.widen(10).reduce { x, y -> x.extra }
        assert wr[0].extra.foo == "hello"
        
        
    }
    
    @Test
    void testIntersect() {
        Region r1 = new Region('chr1', 1, 100)
        Region r2 = new Region('chr1', 1, 120)        
        Region ix = r1.intersect(r2)
        assert ix.size() == 100
    }
    
    @Test
    void testIntersectDifferentChrs() {
        Region r1 = new Region('chr1', 1, 100)
        Region r2 = new Region('chr2', 1, 120)        
        Region ix = r1.intersect(r2)
        assert ix.size() == 0
    }
  
    
    @Test
    void testIntersectNonIntersecting() {
        Region r1 = new Region('chr1', 1, 100)
        Region r2 = new Region('chr1', 150, 220)        
        Region ix = r1.intersect(r2)
        assert ix.size() == 0
    }
    
    @Test
    void 'regions on different chrs have no mutual overlap'() {
        Region r1 = new Region('chr2', 1, 100)
        Region r2 = new Region('chr1', 1, 120)
        approx(r1.mutualOverlap(r2),0d)
     }
     
    @Test
    void 'regions on same chr with no overlap have zero mut overlap'() {
        Region r1 = new Region('chr1', 1, 100)
        Region r2 = new Region('chr1', 200, 220)
        approx(r1.mutualOverlap(r2),0d)
     }
  
    @Test
    void 'regions that are identical have mutual overlap of one'() {
        Region r1 = new Region('chr1', 1, 100)
        Region r2 = new Region('chr1', 1, 100)
        approx(r1.mutualOverlap(r2),1.0d)
     }
     
    @Test
    void 'half overlapping regions have mutual overlap of 50'() {
        Region r1 = new Region('chr1', 1, 100)
        Region r2 = new Region('chr1', 50, 150)
        approx(r1.mutualOverlap(r2),0.5d)
     } 
     
    @Test
    void 'mutual overlap of small region overlapping large'() {
        Region r1 = new Region('chr1', 1, 100)
        Region r2 = new Region('chr1', 48, 50)
        approx(r1.mutualOverlap(r2),0.02d)
     } 
  
      
     // Simplistic but easy to use approximate equals
     void approx(double x1, double x2) {
         assert Math.abs(x2 - x1) < 0.01
     }
}
