package gngs

import static org.junit.Assert.*

class XR extends RegulatingActor<String> {
    
    int waitCount = 0
    
    void process(String m) {
        println "String: $m"
    }
    
    void onEnd() {
        for(i in 0..5) {
            println "Waiting $i"
            waitCount++
            Thread.sleep(1000)
        }
    }
}

import org.junit.Test

class RegulatingActorTest {

    @Test
    public void test() {
        XR xr = new XR()
        xr.start()
        xr.sendTo('hi')
        xr.sendStop()
        
        println "joining"
        xr.join()
        println "Done"
        
        assert xr.waitCount > 2
        
        Thread.sleep(2000)
    }

}
