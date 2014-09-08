import static org.junit.Assert.*;

import org.junit.Test;

class TimeTest {

    @Test
    public void testDump() {
        
        Time.time("foo") {
            Thread.sleep(2000)
        }
        Time.time("bar") {
            Thread.sleep(2000)
        }
        Time.time("foo") {
            Thread.sleep(1000)
        }
        Time.dump()
    }

    @Test
    public void testProgress() {
        ProgressCounter p = new ProgressCounter(withTime:true)
        p.timeInterval = 20
        p.lineInterval = 10
        for(i in 1..100) {
            p.count()
            Thread.sleep(10)
        }
    }
}
