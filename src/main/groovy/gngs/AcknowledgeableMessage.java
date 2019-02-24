package gngs;

import java.util.concurrent.atomic.AtomicInteger;

public class AcknowledgeableMessage {
    
    public AcknowledgeableMessage(Object payload, AtomicInteger counter) {
        this.payload = payload;
        this.acknowledgeCounter = counter;
    }
    
    final public Object payload;
    
    final public AtomicInteger acknowledgeCounter;
}