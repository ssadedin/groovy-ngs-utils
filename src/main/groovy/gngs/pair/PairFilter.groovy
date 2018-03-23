package gngs.pair

import gngs.CompactReadPair
import gngs.SAMRecordPair
import groovy.transform.CompileStatic
import groovyx.gpars.actor.DefaultActor
import htsjdk.samtools.SAMRecord

class PairFilter extends DefaultActor {
    
    String expr
    
    PairFormatter formatter
    
    GroovyShell shell 
    
    Closure filterClosure
    
    PairFilter(PairFormatter formatter, String script) {
        this.formatter = formatter
        shell = new GroovyShell()
        this.filterClosure = shell.evaluate("""{ pair -> 
            $script
       }""")
    }

    @CompileStatic
    void act() {
        
        loop {
            react { Object msg ->
                assert msg != null
                if(msg == "stop") {
                    terminate()
                }
                else {
                    List msgList = (List)msg
                    process((SAMRecordPair)msgList[0], (SAMRecord)msgList[1])
                }
            }
        }
    }
    
    void process(SAMRecordPair pair, SAMRecord r2) {
        Object result = filterClosure(pair)
        if(result == true) {
            formatter << [pair,r2]
        }
    }
}
