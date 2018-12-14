package gngs

import groovy.transform.CompileStatic
import groovy.util.logging.Log
import groovyx.gpars.actor.Actor
import groovyx.gpars.actor.DefaultActor
import groovyx.gpars.actor.impl.MessageStream

import java.util.concurrent.atomic.AtomicInteger

@CompileStatic
class AcknowledgeableMessage {
    
    AcknowledgeableMessage(Object payload, AtomicInteger counter) {
        this.payload = payload
        this.acknowledgeCounter = counter
    }
    
    Object payload
    
    AtomicInteger acknowledgeCounter
}

/**
 * An Actor that exerts back-pressure to regulate the rate of incoming messages.
 * 
 * @author Simon Sadedin
 *
 * @param <T>
 */
@Log
abstract class RegulatingActor<T> extends DefaultActor {
    
    private final int softLimit
    
    private final int hardLimit
    
    final AtomicInteger pendingMessages = new AtomicInteger(0)
    
    /**
     * Count of messages pending in downstream due to this actor
     */
    final AtomicInteger downstreamCounter = new AtomicInteger(0)
    
    public RegulatingActor(RegulatingActor downstream, int softLimit, int hardLimit) {
        this.softLimit = softLimit;
        this.hardLimit = hardLimit;
        this.downstream = downstream;
    }
    

    public RegulatingActor(int softLimit, int hardLimit) {
        super();
        this.softLimit = softLimit;
        this.hardLimit = hardLimit;
    }


    final RegulatingActor downstream
    
    ProgressCounter progress 
    
    @Override
    @CompileStatic
    void act() {
        loop {
            react { msg ->
                if(msg == "stop") {
                    this.onEnd()
                    if(this.progress != null)
                        this.progress.end()
                    log.info "${this.class.name} terminating"
                    this.terminate()
                }
                else {
                    if(this.progress != null)
                        this.progress.count()
                    AcknowledgeableMessage am = (AcknowledgeableMessage)msg
                    process((T)am.payload)
                    this.pendingMessages.decrementAndGet()
                    am.acknowledgeCounter.decrementAndGet()
                }
            }
        }
    }
    
    void onEnd() {
    }
    
    abstract void process(T message)
    
    long throttleWarningMs = 0
    
    @CompileStatic
    void sendDownstream(Object message) {
        this.downstream.send(new AcknowledgeableMessage(message, this.downstreamCounter))
    }
    
    @CompileStatic
    public MessageStream sendTo(Object o) {
        this.send(new AcknowledgeableMessage(o, this.pendingMessages))
    }
    
    @CompileStatic
    public MessageStream send(AcknowledgeableMessage message) {
        
        if(message.acknowledgeCounter.get() > softLimit) {
            long nowMs = System.currentTimeMillis()
            if(nowMs - throttleWarningMs > 30000) {
                log.info "Soft throttling $this due to congestion (pending=$message.acknowledgeCounter/$pendingMessages))"
                throttleWarningMs = nowMs
            }
            Thread.sleep(50)
        }
            
        while(message.acknowledgeCounter.get() > hardLimit) {
            log.info "Hard throttling $this due to congestion (pending=$message.acknowledgeCounter/$pendingMessages))"
            Thread.sleep(3000)
        }
        
        message.acknowledgeCounter.incrementAndGet()
        this.pendingMessages.incrementAndGet()
        super.send(message)
    }
}
