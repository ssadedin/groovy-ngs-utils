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
        testRegions.eachLine {
            String [] parts = it.split(":")
            regions.addRegion(parts[0].trim(), parts[1].split("-")[0].toInteger(), parts[1].split("-")[1].toInteger()+1)
        }
        
        Regions reduced = regions.reduce()
        for(Region r in reduced) {
            assert regions.getOverlaps(r).size() > 0
        }
    }

}
