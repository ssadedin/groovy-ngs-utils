package gngs;

import java.util.concurrent.atomic.AtomicInteger;

public class AcknowledgeableMessage<T> {
    
    public AcknowledgeableMessage(T payload, AtomicInteger counter) {
        this.payload = payload;
        this.acknowledgeCounter = counter;
    }
    
    public T payload;
    
    final public AtomicInteger acknowledgeCounter;
}