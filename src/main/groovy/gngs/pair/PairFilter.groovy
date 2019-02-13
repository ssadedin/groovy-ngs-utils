package gngs.pair

import java.util.List

import gngs.CompactReadPair
import gngs.RegulatingActor
import gngs.SAMRecordPair
import groovy.transform.CompileStatic
import groovyx.gpars.actor.DefaultActor
import htsjdk.samtools.SAMRecord

@CompileStatic
class PairFilter extends RegulatingActor<Paired> {
    
    String expr
    
    GroovyShell shell 
    
    Closure filterClosure
    
    PairFilter(RegulatingActor<Paired> formatter, String script) {
        super(formatter,5000,10000)
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
    public void process(Paired paired) {
        process((SAMRecordPair)paired.r1, paired.r2)
    }
}
