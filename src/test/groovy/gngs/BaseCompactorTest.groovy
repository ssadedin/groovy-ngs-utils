package gngs

import static org.junit.Assert.*

import org.junit.Test

import static gngs.BaseCompactor.*

class BaseCompactorTest {

    @Test
    void 'bases compact to expected bytes'() {
        
        String bases = 'ATCGCGACTCG'
        byte [] compressed = compact(bases.bytes)
        
        println "Compressed: " + compressed
        
        //                       A     T
        assert compressed[0] == (1i | (2i<<4i))
        
        // Because the sequenece has an odd length, the last compacted value should be that base alone
        assert compressed[-1] == 4
    }
    
    @Test
    void 'compressing bases to bytes uncompresses to same bases'() {
        String bases = 'ATCGCGACTCG'
        byte [] compressed = compact(bases.bytes)
        byte [] uncompressed = expand(compressed)
        
        println "Uncompressed: " + uncompressed
        
        assert bases.bytes == uncompressed
        
        assert new String(uncompressed) == bases
        
         
    }
}
