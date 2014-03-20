import static org.junit.Assert.*;

import org.junit.Test;


class StatsTest {

    def x = new Stats()
    
    @Test
    public void test() {
        
        x << 4
        x << 7
        x << 8
        
        println x.mean
        
        println Stats.from(["cat","dog","treehouse","farm"]) { it.size() }
        
        println Stats.from([5,4,2,9])
    }
    
    @Test
    void testMean() {
        
        def values = [1, 5, 7, 8]
        Iterator i = values.iterator()
        
        println Stats.mean { 
           i.next() 
        }
        
        println Stats.mean(values)
        
        println Stats.mean(values.iterator())
    }
    
    @Test
    void testPercentile() {
        
        def values = [1, 5, 7, 8]
        Iterator i = values.iterator()
         println Stats.percentile(200) {
             i.next()
         }.getPercentile(50)
    }

}
