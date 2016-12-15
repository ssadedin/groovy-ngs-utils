import static org.junit.Assert.*;

import org.junit.Before;
import org.junit.Test;

import static RelationshipType.*

class SamplesToPedTest {

    Pedigree p = new Pedigree(id:"testfamily")
        
    Subject mum = new Subject(id: "mum", relationships:[new Relationship(from:"mum", type:MOTHER)])
        
    Subject dad = new Subject(id: "dad", relationships:[new Relationship(from:"dad", type:FATHER)])
        
    Subject daughter = new Subject(id: "daughter", relationships:[new Relationship(from: "daughter", type:DAUGHTER)])
    
    @Before
    void setup() {
        p.individuals=[mum,dad,daughter]
    }
        
    @Test
    public void testLcs() {
        
        SamplesToPed stp = new SamplesToPed()
        
        assert stp.lcs(["fabcd", "fabce","fabcf","frog","dog","hat"]) == [group: "fabc", members:["fabcd","fabce","fabcf"]]
        
        assert stp.lcs(["cow", "tree","house"])== [group:"", members:["cow","tree","house"]]
        
        
        assert stp.lcs(["J00170", "J00170Mic","J00170Mum"]) ==  "J00170"
    }
    
    @Test
    void testLinkParent() {
        SamplesToPed.linkParentRelationship(p, daughter, MOTHER)
        
        assert mum.relationships[0].to == daughter.id
        assert daughter.relationships.find { it.type == DAUGHTER && it.to == mum.id }
        
        SamplesToPed.linkParentRelationship(p, daughter, FATHER)
        
        assert mum.relationships[0].to == daughter.id
        assert dad.relationships[0].to == daughter.id
        assert daughter.relationships.find { it.type == DAUGHTER && it.to == mum.id }
        assert daughter.relationships.find { it.type == DAUGHTER && it.to == dad.id }
    }
    
    @Test 
    void testFillRelationships() {
        SamplesToPed.fillRelationships(p)
        assert mum.relationships[0].to == daughter.id
        assert dad.relationships[0].to == daughter.id
        assert daughter.relationships.find { it.type == DAUGHTER && it.to == mum.id }
        assert daughter.relationships.find { it.type == DAUGHTER && it.to == dad.id }        
    }
}
