import static org.junit.Assert.*

import org.junit.Test

import gngs.PrefixTrie
import gngs.Utils
import trie.TrieNode

class PrefixTrieTest {
    
    PrefixTrie trie = new PrefixTrie()

//    @Test
    public void testAdd() {
        trie.add("ACTG","A")
        
        assert trie.root["A"] != null
        assert trie.root["A"].children["C"] != null
        assert trie.root["A"].children["C"].children["T"] != null
        assert trie.root["A"].children["C"].children["T"].children["G"] != null
        
        assert trie.root["A"].children["C"].children["T"].children["G"].values[0] == "A"
    }
    
//    @Test
    void testGetAt() {
        trie.add("ACTG","cat")
        assert trie.getAt("ACTG") == ["cat"]
        assert trie.getAt("ACT") == []
        assert trie["ACTG"] == ["cat"]
    }
    
    
    @Test
    void testGetByPrefix() {
        trie.add("ACTGTCG","cat")
        
        assert trie.startsWith("ACT") == ["cat"]
        assert trie.startsWith("ACTG") == ["cat"]
        assert trie.startsWith("ATTG") == []
        assert trie.startsWith("ACTGTCG") == ["cat"]
        assert trie.startsWith("ACTGTCT") == []
        assert trie.startsWith("") == ["cat"]
        
        trie.add("ACTCTA","dog")
        assert trie.startsWith("ACT").sort() == ["cat","dog"]
        assert trie.startsWith("ACTG").sort() == ["cat"]
        assert trie.startsWith("ACTC").sort() == ["dog"]
    }
    
    @Test
    void testMismatches() {
        trie.add("ACTGTCG","cat")
        assert trie.startsWith("ACT",1) == ["cat"]
        assert trie.startsWith("TCT",1) == ["cat"]
        assert trie.startsWith("TTT",1) == []
        
        assert trie.startsWith("TTT",2) == ["cat"]
        assert trie.startsWith("TTG",2) == []
    }

    @Test
    void testDeletions() {
        trie.add("ACTGTCG","cat")
        assert trie.startsWith("AT",0,0,1) == ["cat"]
        assert trie.startsWith("ACG",0,0,1) == ["cat"]

    }
    
    @Test
    void testInsertions() {
        trie.add("ACTGTCG","cat")
        assert trie.startsWith("ACCT",0,1,0) == ["cat"]
        assert trie.startsWith("ACTTGTCG",0,1,0) == ["cat"]
    }
    
    @Test
    void testSingleBase() {
        trie.add("A","cat")
        assert trie.startsWith("A") == ["cat"]
    }
    
    @Test
    void testRealSequence() {
        String ref =   "TATTTTTTAATAATTC"    
        String query =      "TTTAATAATTCATATA"
        
        trie.add(ref.reverse(),"cat")
        
        Utils.time("complex Trie query") {
            println trie.startsWith(query.reverse(), 1,1,1,5)
        }
    }
    
    @Test
    void testOrder() {
        trie.add("ACTGTCG","cat")
        trie.add("ACTCTCG","dog")
        trie.add("ATTGTCG","tree")
        
        assert trie.startsWith("ACTCT",2) == ["dog","cat","tree"]
    }
    
    @Test
    void testCost() {
        
       String ref =  "TAAAGTAAAATAGGA"
       String query = "AAAGTAAAATAGGAC"
       
       TrieNode.verbose = true
        
       trie.add(ref, "foo") 
       
       List result = trie.query(query, 1,1,1,1)
       
       println "Result = " + result + " cost = " + result[0].cost(TrieNode.DEFAULT_COSTS)
    }
    
    @Test 
    void testDoubleAdd() {
        trie = new PrefixTrie()
        trie.add("CGGCCTGGGCACGCA","cat")
        trie.add("ACGCACGGGTCCGGC","bar")
        // trie.startsWith("CCCGGACGCACGGG",1,1,1,5)
        assert trie.startsWith("CCCGGACGCACGGG",1,1,1,5).size() == 1
    }
    
    @Test
    void testNodeIterateKeys() {
        trie = new PrefixTrie()
        trie.add("foo",1)
        trie.add("fo",2)
        
        List keys = trie.root["f"].keyIterator("f").collect { it }
        
        println "Keys are: $keys"
        
        assert ["fo","foo"].every { it in keys }
    }
    
    @Test
    void testIterateKeys() {
        trie = new PrefixTrie()
        trie.add("foo",1)
        trie.add("fo",2)
        
        List keys = trie.keyIterator().collect { it }
        
        println "Keys are: $keys"
        
        assert ["fo","foo"].every { it in keys }
    } 
}
