import java.util.concurrent.atomic.AtomicInteger

import groovy.transform.CompileStatic
import groovy.util.logging.Log
import groovyx.gpars.actor.DefaultActor
import htsjdk.samtools.SAMFileWriter

/**
 * An actor that adds asynchronous writing to an OrderedPairWriter
 * <p>
 * Note that if you feed reads directly via the << operator,
 * you may feed them faster than the underlying writer can write, which 
 * will cause the writes to be queued. This queue can expand, taking 
 * a large amount of memory. To limit memory, use the
 * {@link#addAlignmentPair} method which will block if the number of pending
 * writes execeeds {@link#maxPendingWrites}.
 * 
 * @author Simon Sadedin
 */
@Log
class OrderedPairActor extends DefaultActor {
    
    OrderedPairWriter writer
    
    boolean stopped = false
    
    AtomicInteger pending = new AtomicInteger()
    
    int maxPendingWrites = 50000
    
    OrderedPairActor(OrderedPairWriter writer) {
        this.writer = writer
    }
    
    Object lock = new Object()
    
    void addAlignmentPair(SAMRecordPair pair) {
        pending.incrementAndGet()
        while(pending.get()>maxPendingWrites) {
            Thread.sleep(50)
        }
        this << pair
    }

    @CompileStatic
    void act() {
        loop {
            react { msg ->
                if(msg == "stop") {
                    stopped = true
                    synchronized(lock) {
                        lock.notifyAll()
                    }
                    log.info "Actor Stopped" 
                }
                else {
                    writer.addAlignmentPair((SAMRecordPair)msg)
                    pending.decrementAndGet()
                }
            }
        }
    }
    
    void shutdown() {
        this << "stop"
        int n = 0
        while(!stopped)  {
            synchronized(lock) {
                lock.wait(500)
                if(n++ % 100 == 0)
                    log.warning("Actor shutdown has taken > " + (n*0.5) + " seconds")
            }
        }
        log.info "Actor is stopped"
    }
}
