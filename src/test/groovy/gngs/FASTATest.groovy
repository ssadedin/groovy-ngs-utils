package gngs

import static org.junit.Assert.*;

import org.junit.Test;

class FASTATest {

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

    @Test
    void 'bases compact to expected bytes'() {
        
        String bases = 'ATCGCGACTCG'
        byte [] compressed = FASTA.compact(bases.bytes)
        
        println "Compressed: " + compressed
        
        //                       A     T
        assert compressed[0] == (1i | (2i<<4i))
        
        // Because the sequenece has an odd length, the last compacted value should be that base alone
        assert compressed[-1] == 4
    }
    
    @Test
    void 'compressing bases to bytes uncompresses to same bases'() {
        String bases = 'ATCGCGACTCG'
        byte [] compressed = FASTA.compact(bases.bytes)
        byte [] uncompressed = FASTA.expand(compressed)
        
        println "Uncompressed: " + uncompressed
        
        assert bases.bytes == uncompressed
        
        assert new String(uncompressed) == bases
        
         
    }
}
