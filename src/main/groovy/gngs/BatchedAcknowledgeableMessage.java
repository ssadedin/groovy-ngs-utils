package gngs;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

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
public class BatchedAcknowledgeableMessage<T> extends AcknowledgeableMessage {
    
    public List<T> pendingMessages;
   
    private final int batchSize;
    
    public BatchedAcknowledgeableMessage(AtomicInteger counter, final int batchSize) {
        super(new ArrayList<T>(batchSize), counter);
        this.batchSize = batchSize;
        this.pendingMessages = (List<T>) this.payload;
    }
    
    public BatchedAcknowledgeableMessage(final BatchedAcknowledgeableMessage<T> other) {
        super(other.pendingMessages, other.acknowledgeCounter);
        this.batchSize = other.batchSize;
    }
     
    public void batchTo(T o, RegulatingActor<T> to) {
        pendingMessages.add(o);
        if(pendingMessages.size()>batchSize) {
            flush(to);
        }
    }
    
    public void flush(RegulatingActor<T> to) {
        this.acknowledgeCounter.addAndGet(pendingMessages.size()-1);
        to.sendLimited(new BatchedAcknowledgeableMessage<T>(this));
        pendingMessages = new ArrayList<T>();
    }
}

