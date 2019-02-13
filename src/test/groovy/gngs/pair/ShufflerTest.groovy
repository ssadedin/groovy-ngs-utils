package gngs.pair

import static org.junit.Assert.*

import org.junit.Test

import gngs.RegulatingActor
import gngs.SAM
import gngs.SAMRecordPair
import htsjdk.samtools.SAMRecord

class ShufflerTest {

    @Test
    public void test() {
        
        List received = []
        RegulatingActor a = RegulatingActor.actor { msg ->
            println "Message: ${msg[1]}"
            received << msg[1]
        } 
        a.start()
        
        SAM bam = new SAM('src/test/data/small.test.bam')
        
        Shuffler s = new Shuffler(a, 10)
        
        List rawReads = bam.withIterator { it.take(50).collect() }
        Map reads = rawReads.groupBy { it.readName }.grep { it.value.size() == 2 }.collectEntries()
        
        println "Got ${reads.size()} reads"
        
        reads.each { e ->
            SAMRecord r = e.value[0]
            Paired p = new Paired(new SAMRecordPair(r1:r), r)
            s.process(p)
        }
        
        // First read out  should be the smallest of the first 10
        assert received[0].readName == reads*.key[0..10].min()
        
        println reads*.key[0]
        
    }

}
