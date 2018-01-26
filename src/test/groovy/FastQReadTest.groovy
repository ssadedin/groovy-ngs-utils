import static org.junit.Assert.*;

import org.junit.Test;

import gngs.FASTQRead


class FastQReadTest {

    @Test
    public void test() {
        def read = new FASTQRead("@M01158:5:000000000-A5RYT:1:1101:18905:1733 1:N:0:9", "ACGT","AAAA")
        
        assert read.name == "@M01158:5:000000000-A5RYT:1:1101:18905:1733"
    }

}
