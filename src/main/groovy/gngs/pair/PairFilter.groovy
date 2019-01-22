package gngs.pair

import java.util.List

import gngs.CompactReadPair
import gngs.RegulatingActor
import gngs.SAMRecordPair
import groovy.transform.CompileStatic
import groovyx.gpars.actor.DefaultActor
import htsjdk.samtools.SAMRecord

@CompileStatic
class PairFilter extends RegulatingActor<List> {
    
    String expr
    
    PairFormatter formatter
    
    GroovyShell shell 
    
    Closure filterClosure
    
    PairFilter(PairFormatter formatter, String script) {
        super(formatter,5000,10000)
        this.formatter = formatter
        shell = new GroovyShell()
        this.filterClosure = (Closure)shell.evaluate("""{ pair -> 
            $script
       }""")
    }

   
    void process(SAMRecordPair pair, SAMRecord r2) {
        Object result = filterClosure(pair)
        if(result == true) {
            this.sendDownstream([pair,r2])
        }
    }

    @Override
    public void process(List msgList) {
        process((SAMRecordPair)msgList[0], (SAMRecord)msgList[1])
    }
}
