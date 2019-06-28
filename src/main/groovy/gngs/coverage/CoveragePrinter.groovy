package gngs.coverage

import gngs.ProgressCounter
import gngs.RegulatingActor
import groovy.transform.CompileStatic
import groovy.util.logging.Log

import java.text.NumberFormat
import java.util.concurrent.ConcurrentLinkedQueue

@CompileStatic
@Log
class CoveragePrinter extends Thread {
    
    Writer w
    
    StringBuilder line = new StringBuilder(2048)
    
    ProgressCounter progress = new ProgressCounter(timeInterval: 10000, lineInterval: 500, log:log)
    
    ConcurrentLinkedQueue<CoverageSummary> queue = new ConcurrentLinkedQueue()
    
    final NumberFormat numberFormat = NumberFormat.numberInstance
    
    CoveragePrinter(Writer w) {
        this.w = w
        this.numberFormat.groupingUsed = false
    }
    
//    @Override
    public void process(final CoverageSummary summary) {
        writePosition(summary)
        progress.count()
    }
    
    boolean continueRunning = true
    
    @Override
    public void run() {
        log.info "Starting printer queue"
        while(continueRunning) {
            CoverageSummary summary = queue.poll()
            if(!summary.is(null)) {
                process(summary)
            }
            else {
                Thread.sleep(100)
            }
        }
    }

    @CompileStatic
    private void writePosition(final CoverageSummary summary) {
        
       if(w == null)
           return
        
       final PositionCounts countInfo = summary.countInfo
       final List<Double> values = summary.values
       final Double coeffV  = summary.coeffV
       
        line.setLength(0)
            
        w.write(countInfo.chr)
        w.write('\t')
        w.write(countInfo.pos.toString())
        w.write('\t')
        if(coeffV != null) {
            w.write(numberFormat.format(coeffV))
            w.write('\t')
        }
        final int numValues = values.size()-1
        int i = 0;
        for(;i<numValues; ++i) {
            w.write(numberFormat.format(values[i]))
            w.write('\t')                                
        }
        w.write(numberFormat.format(values[numValues]))
        w.write('\n')
    }
    
}
