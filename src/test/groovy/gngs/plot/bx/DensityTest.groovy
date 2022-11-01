package gngs.plot.bx

import static org.junit.Assert.*

import java.awt.Color

import org.junit.Test

class DensityTest {

    Random r = new Random()

    @Test
    public void test() {
        def x = (1..1000).collect { r.nextGaussian() }
        def pd = new Density.Points(data:x, color: Color.green)
        
        println "Done"
    }
    
    @Test
    void 'Test all zeros'() {
        def x = (1..10).collect {0.0d}
        def pd = new Density.Line(data:x, color: Color.green)
         
        println "OK"
    }
    
    @Test
    void 'test from double array'() {
        double[] x =  (1..1000).collect { r.nextGaussian() } as double[]
        def pd = new Density.Line(data:x, color: Color.green)
        
        assert pd.step > 0.0d
        assert pd.x.size() > 10
    }

}
