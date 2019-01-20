package gngs;

import static org.junit.Assert.*

import org.junit.Test

class CoverageGapsTest {

    @Test
    public void testSimpleGap() {
        
        Utils.configureSimpleLogging()
        
        List covs = createCov()
        
        // Create a region of low coverage from 130 - 134
        covs[30][4] = 19
        covs[31][4] = 18
        covs[32][4] = 18
        covs[33][4] = 19
        
        CoverageGaps gaps = calc(covs)
        
        assert gaps.blocks.size() == 1
        
        gaps.blocks[0].with {
            assert start == 130
            assert end == 133
            assert size() == 4
        }
    }
    
    @Test
    void testAdjacentGaps() {
        List covs = createCov()
        
        // ----____-___-----
        (30..33).each { covs[it][4] = 19 }
        (35..37).each { covs[it][4] = 19 }
        
        CoverageGaps gaps = calc(covs)
        assert gaps.blocks.size() == 2
        
        gaps.with {
            assert blocks[0].start == 130
            assert blocks[0].end == 133
            
            assert blocks[1].start == 135
            assert blocks[1].end == 137
            assert blocks[1].size() == 3
        }
    }
    
    @Test
    void test1BaseGap() {
        List covs = createCov()
        
        // ----_------
        covs[30][4] = 19
        
        CoverageGaps gaps = calc(covs)
        
        assert gaps.blocks.size() == 1
        
        gaps.with {
            assert blocks.size() == 1
            assert blocks[0].start == 130
            assert blocks[0].end == 130
            assert blocks[0].size() == 1
        }
    }
    
    List createCov() {
        List covs = []
        for(i in 1..100)
            covs << ["chr1",100,200, 0, 25] 
        covs.eachWithIndex { c,i -> c[3] = i }
        return covs
    }
    
    CoverageGaps calc(List covs) {
        CoverageGaps gaps = new CoverageGaps("test.cov")
        new ByteArrayInputStream(covs*.join('\t').join('\n').bytes).withReader { r ->
            gaps.calculateFromBEDTools(r)
        }
        return gaps
    }
}
