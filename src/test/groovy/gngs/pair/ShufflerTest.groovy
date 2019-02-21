package gngs.pair

import static org.junit.Assert.*

import org.junit.Test

import gngs.RegulatingActor
import gngs.SAM
import gngs.SAMRecordPair
import htsjdk.samtools.SAMRecord

class ShufflerTest {

    @Test
    public void 'test randomised order'() {
        
        List received = []
        RegulatingActor a = RegulatingActor.actor { msg ->
            println "Message: ${msg.r1}"
            received << msg.r1
        } 
        a.start()
        
        SAM bam = new SAM('src/test/data/small.test.bam')
        
        Shuffler s = new Shuffler(a, 10)
        
        List rawReads = bam.withIterator { it.take(50).collect() }
        Map reads = rawReads.groupBy { it.readName }.grep { it.value.size() == 2 }.collectEntries()
        
        println "Got ${reads.size()} reads"
        
        reads.each { e ->
            SAMRecord r = e.value[0]
            Paired p = new Paired(r.readName, new SAMRecordPair(r1:r), new SAMRecordPair(r1:r), false, true, 0)
            s.process(p)
        }
        
        // the order should be randomised!
        assert received[0..10]*.readName != rawReads[0..10]*.readName
        
        println reads*.key[0]
        
    }

}
