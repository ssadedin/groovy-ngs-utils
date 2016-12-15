import static org.junit.Assert.*;

import org.junit.Test;


class OMIMTest {

    @Test
    public void testParse() {
        
        def omim = OMIM.parse("/Users/simon/Downloads/omim/genemap2.txt")
        
        assert omim.size() > 10
        
        assert omim["ANIB3"].disorders.any { it.contains("Aneurysm") }
        
//        assert omim["INPP5D"].disorders.any { it.contains("Jouber") }
        
        
    }

}
