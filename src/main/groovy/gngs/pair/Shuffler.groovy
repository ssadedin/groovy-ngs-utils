package gngs.pair

import java.util.concurrent.ConcurrentSkipListMap
import java.util.concurrent.PriorityBlockingQueue
import java.util.concurrent.atomic.AtomicInteger
import gngs.ProgressCounter
import gngs.RegulatingActor
import groovy.transform.CompileStatic
import groovy.util.logging.Log
import groovyx.gpars.actor.impl.MessageStream

@CompileStatic
class PairedReadComparator implements Comparator<Paired> {
    @Override
    public int compare(Paired p1, Paired p2) {
        return p1.flags.compareTo(p2.flags)
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
    
    int emitThreshold
    
    static NavigableMap buffer  = new ConcurrentSkipListMap<Integer, Paired>()
    
    static int shufflerCounter = 0
    
    Random random = new Random(shufflerCounter++)
    
    boolean threadStarted = false
    
    int minDistance = 5000
    
    static AtomicInteger currentBufferSize = new AtomicInteger(0)
    
    Shuffler(RegulatingActor downstream, int bufferSize) {
        super(downstream, 5000, (int)10000)
//        this.buffer = new PriorityBlockingQueue(bufferSize, new PairedReadComparator())
        this.bufferSize = bufferSize
        this.emitThreshold = (int)(bufferSize * 0.7d)
        log.info "Shuffler emitting at ${emitThreshold}/${bufferSize} messages"
    }
    
    @CompileStatic
    class EmitThread implements Runnable {
        @Override
        void run() {
            
            ProgressCounter stopProgress = new ProgressCounter(extra: {
                "Writing residual reads from shuffler buffer (${currentBufferSize.get()} remaining)"
            })
            
            while(true) {
                while(Shuffler.this.currentBufferSize.get() > emitThreshold) {
                    emit()
                    if(stopped) {
                        stopProgress.count()
                    }
                }
                
                if(stopped) {
                    stopProgress.end()
                    break
                }
                Thread.sleep(5)
            }
            log.info "Exiting shuffler emit thread"
        }
    }
    
    @Override
    public void process(Paired message) {
        
        if(stopped)
            throw new IllegalStateException("I was stopped!")
        
        if(!threadStarted) {
            log.info "Starting shuffler background thread"
            new Thread(new EmitThread()).start()
            threadStarted = true
        }
        
        int key = random.nextInt()
        while(buffer.containsKey(key)) {
            key = random.nextInt()
        }
        buffer.put(key, message)
        currentBufferSize.incrementAndGet()
        while(currentBufferSize.get() > bufferSize) {
            Thread.sleep(5)
        }
    }
    
    void emit() {
        int index = random.nextInt()
        Map.Entry<Integer,Paired> entry = this.buffer.floorEntry(index)
        if(entry.is(null)) {
            entry = this.buffer.ceilingEntry(index)
        }
        
        // This can happen because the buffer is shared between many threads
        if(entry.is(null)) {
            return
        }

        Paired msg = this.buffer.remove(entry.key)
        if(!msg.is(null)) { // there is a miniscule chance the other thread grabbed the same read pair
            this.currentBufferSize.decrementAndGet()
            this.sendDownstream(msg)
        }
    }

    @Override
    public void onEnd() {
        this.emitThreshold = 0
        
        if(!threadStarted) {
            // Only wait for thread to flush reads if we actually started a thread
            return
        }
            
        while(this.currentBufferSize.get()>0) {
            Thread.sleep(50)
        }
    }
}
