package gngs.coverage

import gngs.*
import groovyx.gpars.actor.Actor
import groovyx.gpars.actor.DefaultActor

import static org.junit.Assert.*

import org.junit.Test

class CoverageCalculatorTest {
    
    static {
        Utils.configureSimpleLogging()
    }


    @Test
    public void test() {
        
        Regions regions = new Regions([
            new Region("chr1:1,338,100-1,338,133")
        ])
        
        List covs = []

        RegulatingActor a = RegulatingActor.actor { SampleReadCount src ->
            println(src)
            covs << src
        }
        
        a.start()
        
        SAM bam = new SAM('src/test/data/test.low.cov.bam')
        CoverageCalculatorActor.processBAM(bam, regions, a, 0, bam.samples[0], false)
        
        a.sendStop()
        a.join()
        
        assert covs.find { it.pos == 1338114 }.reads == 12
        assert covs.find { it.pos == 1338117 }.reads == 12
        assert covs.find { it.pos == 1338118 }.reads == 13
    }
}
