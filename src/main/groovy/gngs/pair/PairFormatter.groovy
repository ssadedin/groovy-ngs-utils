package gngs.pair

import gngs.CompactReadPair
import gngs.SAMRecordPair
import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.util.logging.Log
import groovyx.gpars.actor.DefaultActor
import groovyx.gpars.util.DefaultMessageQueue
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SAMUtils
import htsjdk.samtools.util.SequenceUtil

@CompileStatic
@Log
class PairFormatter extends DefaultActor {
    
    StringBuilder buffer
    
    int maxBufferSize
    
    PairWriter writer
    
    int formatted = 0
    
    int buffered = 0
    
    /**
     * If true, the original position of the read will be appended 
     * to the read name. This allows backtracking of reads to compare between
     * alignments, etc.
     */
    boolean addPosition = false
    
    PairFormatter(int bufferSize, PairWriter writer) {
        buffer = new StringBuilder(bufferSize+2000)
        this.maxBufferSize = bufferSize
        assert writer != null
        this.writer = writer
    }

    @CompileDynamic
    void act() {
        
        loop {
            react { Object msg ->
                assert msg != null
                if(msg == "stop") {
                    terminate()
                }
                else {
                    List msgList = (List)msg
                    process((CompactReadPair)msgList[0], (SAMRecord)msgList[1])
                }
                
            }
        }
    }
    
    void process(CompactReadPair pair, SAMRecord r2) {
        pair.appendTo(buffer, r2, this.addPosition)
        assert writer != null
        assert buffer != null
        
        if(buffer.size() > maxBufferSize) {
            writer << [ content: buffer.toString(), reads: buffered ]
            writer.pending.addAndGet(buffered)
            buffer.setLength(0)
            formatted += buffered
            buffered = 0
        }
        else {
            buffered += 2
        }
    }
}
