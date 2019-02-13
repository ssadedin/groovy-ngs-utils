package gngs.pair

import java.util.concurrent.PriorityBlockingQueue

import gngs.RegulatingActor
import groovy.transform.CompileStatic
import groovy.util.logging.Log
import groovyx.gpars.actor.impl.MessageStream

@CompileStatic
class ReadNameComparator implements Comparator<Paired> {
    @Override
    public int compare(Paired p1, Paired p2) {
        return p1.r2.readName.compareTo(p2.r2.readName)
    }
}

/**
 * An actor which stores a buffer of its messages and emits the incoming messages in random order
 * 
 * @author Simon Sadedin
 */
@CompileStatic
@Log
class Shuffler extends RegulatingActor<Paired> implements Runnable {
    
    final int bufferSize
    
    final int emitThreshold
    
    PriorityBlockingQueue buffer 
    
    Random random = new Random(0)
    
    boolean threadStarted = false
    
    int minDistance = 5000
    
    Shuffler(RegulatingActor downstream, int bufferSize) {
        super(downstream, 5000, (int)10000)
        this.buffer = new PriorityBlockingQueue(bufferSize, new ReadNameComparator())
        this.bufferSize = bufferSize
        this.emitThreshold = (int)(bufferSize * 0.7d)
        log.info "Shuffler emitting at ${emitThreshold}/${bufferSize} messages"
    }
    
    @CompileStatic
    class EmitThread implements Runnable {
        @Override
        void run() {
            while(true) {
                while(Shuffler.this.buffer.size() > emitThreshold) {
                    emit()
                    if(stopped) {
                        log.info "Shuffler stopping emit thread due to stop flag set"
                    }
                }
                Thread.sleep(5)
                
                if(stopped)
                    break
            }
        }
    }

    @Override
    public void process(Paired message) {
        
        if(!threadStarted) {
            log.info "Starting shuffler background thread"
            new Thread(new EmitThread()).start()
            threadStarted = true
        }
        
        buffer.add(message)
        while(buffer.size() > bufferSize) {
            Thread.sleep(5)
        }
    }
    
    void emit() {
        int index = (int)(Math.round(random.nextDouble() * (this.buffer.size()-1)))
        Paired msg = this.buffer.take()
        this.sendDownstream([msg.r1, msg.r2])
    }

    @Override
    public void onEnd() {
        while(!this.buffer.isEmpty())
            this.emit()
    }
}
