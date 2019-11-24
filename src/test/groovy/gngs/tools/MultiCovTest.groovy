package gngs.tools

import static org.junit.Assert.*

import graxxia.Matrix
import graxxia.TSV
import org.junit.Test

class MultiCovTest {
    
    MultiCov mc = new MultiCov()
    
    String BED = 'src/test/data/multicov_test.bed'
    
    String BAM = 'src/test/data/small.test.bam'

    @Test
    public void 'Test basic coverage calculation'() {
        def outputFile= 'src/test/data/small.test.bam.cov.txt'
        mc.test(['-bed', BED, '-o',outputFile,BAM]) {
            mc.run()
        }
        
        Matrix covs = new Matrix(new TSV(outputFile, columnNames:['chr','pos','cov']))
        assert covs.grep { pos == 2682 }.cov[0] == 25
        assert covs.grep { pos == 2772 }.cov[0] == 7
    }
    
    @Test
    void 'Same sample in multiple BAMs'() {
        def outputFile= 'src/test/data/small.test.bam.cov.txt'
        mc.test(['-bed', BED, '-o',outputFile,BAM,BAM]) {
            mc.run()
        }
        
        Matrix covs = new Matrix(new TSV(outputFile, columnNames:['chr','pos','cov1', 'cov2']))
        
        assert covs.grep { pos == 2682 }.cov1[0] == 25
        assert covs.grep { pos == 2682 }.cov2[0] == 25
        assert covs.grep { pos == 2772 }.cov1[0] == 7
        assert covs.grep { pos == 2772 }.cov2[0] == 7
    }
}
