package gngs.pair

import java.io.Writer
import java.util.concurrent.atomic.AtomicInteger
import gngs.SAMRecordPair
import groovy.transform.CompileStatic
import groovyx.gpars.actor.DefaultActor

@CompileStatic
class PairWriter extends DefaultActor {
    
    Writer out
    
    int written = 0
    
    public AtomicInteger pending = new AtomicInteger()
    
    PairWriter(Writer writer) {
        this.out = writer
    }
    
    void act() {
        loop {
            react { msg ->
                if(msg == "stop")
                    terminate()
                else {
                    Map msgMap = (Map)msg
                    out.append((String)msg['content'])
                    int readCount = (int)msg['reads']
                    written = written + readCount
                    pending.addAndGet(-readCount)
                }
            }
        }
    }
}


