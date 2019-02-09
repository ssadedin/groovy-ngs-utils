package gngs.tools

import static org.junit.Assert.*

import htsjdk.variant.variantcontext.VariantContext
import org.junit.Test

class CreatePopulationStatisticsVCFTest {

    @Test
    public void test() {
        
        List<Map> results = []
        CreatePopulationStatisticsVCF cpsv = new CreatePopulationStatisticsVCF() {
            void printVCFSite(VariantContext v0, int ac, int an, int gtc) {
                results << [vc: v0, ac: ac, an: an, gtc: gtc]
            }
        }
        
        cpsv.test(['-o', '/dev/null','src/test/data/trio.test.vcf.gz']) { 
            
            cpsv.run()
            
            assert cpsv.sampleSexes['NA12878'] == gngs.Sex.FEMALE
            assert cpsv.sampleSexes['NA12877'] == gngs.Sex.MALE
            assert cpsv.sampleSexes['NA12879'] == gngs.Sex.FEMALE
            
            println "Calculated ${results.size()} results"
            
            assert results.every { it.gtc == 3 }
            assert results*.an.take(100).every { it == 5 }
            
            Map v1 = results.find { it.vc.start == 3241317 }
            assert v1.ac == 3
            
            Map v2 = results.find { it.vc.start == 2853146 }
            assert v2.ac == 2 
            
            Map v3 = results.find { it.vc.start == 2825403 }
            assert v3.ac == 1
            
            Map v4 = results.find { it.vc.start == 21154426 } // chrY
            assert v4.ac == 1
            assert v4.an == 1
            assert v4.gtc == 3
             
        }
    }

}
