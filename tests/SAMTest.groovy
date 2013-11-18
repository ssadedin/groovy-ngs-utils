import static org.junit.Assert.*;

import org.junit.Test;


class SAMTest {

    @Test
    public void testHardClipped() {
        
        SAM sam = new SAM("tests/data/hard_clipped.bam")
        
        sam.pileup("chr12", 53819703, 53819708) { PileupIterator.Pileup p ->
            println p.position + " : " + p.alignments.size()
        }
        
//        sam.coverage("chr9",98211125,98211141) {  p ->
//             println p.position   
//        }
//        
//        
//        println sam.meanCoverage("chr9", 98211125, 98211141)
        
        println sam.meanCoverage("chr12", 53819703, 53819708)
    }

} 
