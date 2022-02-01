package gngs

import static org.junit.Assert.*

import org.junit.Test

class FastQTest {

    @Test
    public void test() {
        
        List names = []
        List seps = []
        
        FASTQ.eachPair('src/test/data/sra_fastq_r1.fastq.gz', 'src/test/data/sra_fastq_r1.fastq.gz') { r1, r2 ->
            
            assert r1.name == r2.name
            
            names << r1.name
            seps << r1.sep
        }
        
        assert names.size() == 5
        assert seps.size() == 5
        
        assert !seps.any { it == '+' }
    }

}
