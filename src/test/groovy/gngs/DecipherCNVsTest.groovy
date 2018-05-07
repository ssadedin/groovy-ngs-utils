package gngs

import static org.junit.Assert.*

import org.junit.Test

class DecipherCNVsTest {
    
    DecipherCNVs ddd = new DecipherCNVs('src/test/data/ddd_test.txt.gz').parse()

    @Test
    public void testOverlaps() {
        List  cnvs = ddd.queryOverlapping(new Region("chr1:10000-15000"))
        assert cnvs.size()  == 1
        
        cnvs = ddd.queryOverlapping(new Region("chr1:10000-50000"))
        assert cnvs.size()  == 2
    }
    
    @Test
    public void testMaxFreq() {    
        approx ddd.maxFreq(new Region("chr1:12000-14000")),  0.20
        approx ddd.maxFreq(new Region("chr1:12000-54000")),  0.675
    }
    
    void approx(double value, double expected) {
        assert Math.abs(value - expected) < 0.01 : "Value $value is not approximately equal to $expected (difference=${Math.abs(value - expected)})"
    }
}
