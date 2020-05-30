package gngs

import static org.junit.Assert.*;

import org.junit.Test;

class FASTATest {
    
    
    @Test
    void 'test repeat detection'() {
         def r = FASTA.repeatAt('ATATATAT'.bytes, 2)
         println r
         assert r.repetitions == 4
         assert r.motif == 'AT'.bytes
         assert r.position == 0
         
         r = FASTA.repeatAt('ATCTATAT'.bytes, 2)
         println r
         assert r.repetitions == 2
         assert r.motif == 'AT'.bytes
         assert r.position == 4
         
         r = FASTA.repeatAt('ATCATCTAT'.bytes, 3)
         println r
         assert r.repetitions == 2
         assert r.motif == 'ATC'.bytes
         
         // When there is no repeat, should return null
         r = FASTA.repeatAt('ATCACCGAT'.bytes, 3)
         assert r == null
         
         r = FASTA.repeatAt('CATATATAT'.bytes, 2, 1)
         println r
         assert r.repetitions == 4
         assert r.motif == 'AT'.bytes
         
         r = FASTA.repeatAt('AAAAA'.bytes, 1)
         println r
         assert r.repetitions == 5
         assert r.motif == 'A'.bytes
         
         r = FASTA.repeatAt('GTGCTGTAGCT'.bytes, 2)
         println r
         assert r == null

     }

////    @Test
//    public void testHeaders() {
//        FASTA fa = new FASTA("/Users/simon/work/dsd/batch1/mytrim/david/test.fasta")
//        
//        println fa.basesAt("1", 10000, 20000)
//    }
//    
//    @Test
//    public void testBaseCounts() {
//        FASTA fa = new FASTA("/Users/simon.sadedin/work/cpipe/hg19.decoy/ucsc.hg19.with_decoy.fasta")
//        Utils.time('count') {
//            println fa.baseCounts("chr1", 10000, 20000)
//        }
//        
//        Utils.time('gc') {
//            println fa.gc("chr1", 10000, 20000)
//        } 
//    }

}
