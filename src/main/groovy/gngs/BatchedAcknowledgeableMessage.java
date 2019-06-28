package gngs;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

import groovy.transform.CompileStatic;

/**
 * A batched message automatically accumulates a fixed number of messagse before sending 
 * on to a {@link RegulatingActor}.
 * <p>
 * To use a BatchedAcknowledgeableMessage, create it up front, then call the 
 * {@link #batchTo(Object, RegulatingActor)} method to send messages. This will
 * accumulate up to {@link #batchSize} messages before actually placing in the queue.
 * <p>
 * When sufficient messages are accumulated, they are all sent together.
 * 
 * @author simon.sadedin
 *
 * @param <T>
 */
@CompileStatic
public class BatchedAcknowledgeableMessage<T> extends AcknowledgeableMessage<List<T>> {
    
    private final int batchSize;
    
    public BatchedAcknowledgeableMessage(AtomicInteger counter, final int batchSize) {
        super(new ArrayList<T>(batchSize), counter);
        this.batchSize = batchSize;
    }
    
    public BatchedAcknowledgeableMessage(final BatchedAcknowledgeableMessage<T> other) {
        super(other.payload, other.acknowledgeCounter);
        this.batchSize = other.batchSize;
    }
     
    public void batchTo(T o, RegulatingActor<T> to) {
        payload.add(o);
        if(payload.size()>batchSize) {
            flush(to);
        }
    }
    
    public void flush(RegulatingActor<T> to) {
        this.acknowledgeCounter.addAndGet(payload.size()-1);
        to.sendLimited(new BatchedAcknowledgeableMessage<T>(this));
        this.payload = new ArrayList<T>(batchSize);
    }
}

