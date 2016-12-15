import static org.junit.Assert.*;

import org.junit.Test;

class RegressionTest {

    @Test
    public void testEquation() {
       Regression r = new Regression() 
       r.model {
           y ~ x + z
       }
       
       def x = new Matrix([[2,4,6,8,10]]).transpose()
       def z = new Matrix([[1,8,10,2,5]]).transpose()
       
       println "X = $x"
       
       // Generate our data directly from the model
       def y = 2 * x + z
       println "Y = $y"
       
       Matrix X = new Matrix(x[][0],z[][0])
       X.names = ["x","z"]
       
       r.run(y[][0], X)
       
       println "Model is " + r
       
       println r.estimateRegressionParameters()
    }
}
