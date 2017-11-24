import static org.junit.Assert.*;

import org.junit.Test;

import gngs.*

class FASTATest {

//    @Test
    public void testHeaders() {
        FASTA fa = new FASTA("/Users/simon/work/dsd/batch1/mytrim/david/test.fasta")
        
        println fa.basesAt("1", 10000, 20000)
    }
    
    @Test
    public void testBaseCounts() {
        FASTA fa = new FASTA("/Users/simon.sadedin/work/cpipe/hg19.decoy/ucsc.hg19.with_decoy.fasta")
        Utils.time('count') {
            println fa.baseCounts("chr1", 10000, 20000)
        }
        
        Utils.time('gc') {
            println fa.gc("chr1", 10000, 20000)
        } 
    }

}
