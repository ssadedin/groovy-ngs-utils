import static org.junit.Assert.*;

import org.junit.Test;

class FASTATest {

    @Test
    public void testHeaders() {
        FASTA fa = new FASTA("/Users/simon/work/dsd/batch1/mytrim/david/test.fasta")
        
        println fa.basesAt("1", 10000, 20000)
    }

}
