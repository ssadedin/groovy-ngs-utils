import static org.junit.Assert.*;

import org.junit.Test;


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
    }
    
    @Test
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
}
