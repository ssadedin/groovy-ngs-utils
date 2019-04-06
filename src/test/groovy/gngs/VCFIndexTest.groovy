package gngs

import static org.junit.Assert.*

import org.junit.Test

class VCFIndexTest {

    @Test
    public void 'test that variants are returned for queries overlapping breakpoints'() {
        
        VCFIndex vcf = new VCFIndex('src/test/data/sv.vcf.gz')
        
        
        // chr1 1313940 INV00000105 C   <INV>   .   LowQual IMPRECISE;SVTYPE=INV;SVMETHOD=EMBL.DELLYv0.8.1;CHR2=chr1;END=1315090
        //
        // .... 1313940 <---1kb----> 1315090
        
        int end = 1315090
        int start = 1313940
        
        // Spanning start breakpoint
        List results = vcf.iterator('chr1', start-20, start+20).collect { it }
        println results
        assert results.size() > 0
        assert results.any { it.pos == start }
        
        // Spanning end breakpoint
        results = vcf.iterator('chr1', end-100, end+100).collect { it }
        assert results.size() > 0
        assert results.any { it.pos == start }
        println results
        
        // Overlapping both breakpoints
        results = vcf.iterator('chr1', start-100, end+100).collect { it }
        assert results.size() > 0
        assert results.any { it.pos == start }
        println results
        
        // Overlapping middle, but no breakpoints
        results = vcf.iterator('chr1', start+200, start+300).collect { it }
        assert results.size() > 0
        assert results.any { it.pos == start }
        println results
        
    }

}
