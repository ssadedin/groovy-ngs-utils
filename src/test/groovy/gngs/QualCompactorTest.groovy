package gngs

import static org.junit.Assert.*

import org.junit.Test

class QualCompactorTest {

    @Test
    public void testCompact() {
        String quals = "AAFFFJJJJJJJJJJJJJJJJJJJJJ"
        
        byte [] compacted = QualCompactor.compact(quals.bytes)
        List qs = quals as List
        assert compacted.length == 3 * 2 + 1
        
        quals = "AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFFJJFJJJJJJJJJJJJJJAFFJJJJJJJJJJJJJJJFJJJJ"
        
        println "Size is ${quals.size()}"
        
        compacted = QualCompactor.compact(quals.bytes)
        qs = quals as List
        assert compacted.length == 12 * 2 + 1
  
    }

    @Test
    public void testExpand() {
        String quals = "AAFFFJJJJJJJJJJJJJJJJJJJJJ"
        byte [] compacted = QualCompactor.compact(quals.bytes)
        byte [] expanded = QualCompactor.expand(compacted)
        
        println "original: " + quals.bytes
        println "expanded: " + expanded
        
        assert quals.bytes == expanded
        
    }
}
