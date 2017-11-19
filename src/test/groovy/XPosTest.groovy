import org.junit.Test

import gngs.LXPos

class XPosTest {
    
    @Test
    void testLXpos() {
        
        
        assert LXPos.computePos('chr1',10) ==           10
        assert LXPos.computePos('chr10',10) ==  1000000010
        assert LXPos.computePos('chr11',10) ==  2000000010
        assert LXPos.computePos('chr19',10) == 10000000010
        assert LXPos.computePos('chr2',10) ==  11000000010
        assert LXPos.computePos('chr20',10) ==  12000000010
        assert LXPos.computePos('chr21',10) ==  13000000010
        assert LXPos.computePos('chr3',10) ==   15000000010
        assert LXPos.computePos('chr9',10) ==   21000000010
        assert LXPos.computePos('chrY',10) ==   23000000010
        
        
    }

}
