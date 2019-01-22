package gngs

import groovy.transform.CompileStatic
import groovy.transform.stc.ClosureParams
import groovy.transform.stc.FirstParam
import groovy.transform.stc.FromAbstractTypeMethods
import groovy.transform.stc.SimpleType
import groovy.util.logging.Log
import groovyx.gpars.actor.Actor
import groovyx.gpars.actor.Actors
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
 * <p>
 * The regulating actor counts how many messages are pending vs how many it has processed.
 * If the number of pending messages exceeds a soft threshold, it blocks for a short time
 * before queueing the next message. If the number of messages exceeds a hard threshold,
 * it blocks for a long time before queuing the message. This allows effective control over how
 * much memory the actor is using in its queue.
 * <p>
 * The RegulatingActor tracks how many messages are pending by incremenging a counter
 * ({@link #pendingMessages}). A user of a RegulatingActor should therefore send messages to it 
 * using the {@link #sendTo} method, so that this counter can be incremented.
 * <p>
 * Sometimes you may want a 'per-client' pending count; that is you don't want one over-zealous
 * producer to block other producers. This can be especially important if there are dependencies 
 * such that blocking some producers could result in a deadlock. This can be implemented
 * by supplying your own pendingCount instead of using the shared one. To do that, 
 * you can use the {@link #send(AcknowledgeableMessage message)} method, where you construct the
 * {@link AcknowledgeableMessage} yourself and supply it with your own counter.
 * <p>
 * Since chaining together regulated actors is a common pattern, this class provides support for a
 * default "downstream". If the downstream actor is supplied in the constructor, you can send 
 * messages to it that are automatically per-client limited using the {@link #sendDownstream} 
 * method.
 * 
 * @author Simon Sadedin
 *
 * @param <T> the type of the messages to be processed
 */
@Log
abstract class RegulatingActor<T> extends DefaultActor {
    
    private final int softLimit
    
    private final int hardLimit
    
    /**
     * Count of messages that have been sent to this actor but not processed
     */
    final AtomicInteger pendingMessages = new AtomicInteger(0)
    
    boolean stopped = false
    
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

    /**
     * Create a RegulatingActor with basic defaults for the queue sizes that suit many tasks
     */
    public RegulatingActor() {
        this(1000,5000)
    }

    final RegulatingActor downstream
    
    ProgressCounter progress 
    
    final static Object STOP = new Object()
    
    @Override
    @CompileStatic
    void act() {
        loop {
            react { msg ->
                if(msg.is(STOP)) {
                    this.stopped = true
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
    
    @CompileStatic
    void sendStop() {
        this << RegulatingActor.STOP
    }
    
    abstract void process(T message)
    
    long throttleWarningMs = 0
    
    @CompileStatic
    void sendDownstream(Object message) {
        if(!this.downstream.is(null))
            this.downstream.send(new AcknowledgeableMessage(message, this.downstreamCounter))
    }
    
    @CompileStatic
    public MessageStream sendTo(T o) {
        this.send(new AcknowledgeableMessage(o, this.pendingMessages))
    }
    
    @CompileStatic
    public MessageStream send(AcknowledgeableMessage message) {
        sendLimited(message)
    }
    
    @CompileStatic
    public MessageStream sendLimited(AcknowledgeableMessage message) {
        
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
    
    @CompileStatic
    static <T> RegulatingActor<T> actor(@ClosureParams(value=FromAbstractTypeMethods) Closure c) {
        RegulatingActor ds = new RegulatingActor<T>() {
            void process(T value) {
                c(value)
            }
        }
        return ds
    }
}
