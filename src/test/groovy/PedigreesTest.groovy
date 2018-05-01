import static org.junit.Assert.*;

import org.junit.Test;

import gngs.Pedigrees


class PedigreesTest {

    @Test
    void testGetFemales() {
        Pedigrees peds = Pedigrees.parse("src/test/data/test.ped")
        
        assert peds.females.sort() == [
            "IND_FAM1_3",
            "IND_FAM2_4",
            "IND_FAM2_6",
            "IND_FAM3_9",
            "IND_FAM4_12",
            "IND_FAM5_16",
            "IND_FAM6_19",
            "IND_FAM7_22"
        ]
    }
    
    @Test
    void testParentsFirst() {
        List pedData = [
            ["F","A",0,0,1,1],
            ["F","B",0,0,2,1],
            ["F","C","A","B",1,1],
        ]
        
        String pedString = pedData*.join('\t').join('\n')+'\n'

        println "PED file: \n" + pedString
        
        def ped = Pedigrees.parse(new StringReader(pedString))
        
        assert ped.motherOf("C").id == "B"
    }
    
    @Test
    void testFailParentWrongSex() {
        List pedData = [
            ["F","A",0,0,1,1],
            ["F","B",0,0,1,1],
            ["F","C","A","B",1,1],
        ]
        
        String pedString = pedData*.join('\t').join('\n')+'\n'

        println "PED file: \n" + pedString
        
        boolean caught = false
        try {
            def ped = Pedigrees.parse(new StringReader(pedString))
        }
        catch(IllegalStateException e) {
            // expected
            caught = true
        }
        
        assert caught : "Failed to throw exception when mother was specified as male"
    } 
    
    @Test
    void testRenameFamily() {
        List pedData = [
            ["F","A",0,0,1,1],
            ["F","B",0,0,2,1],
            ["F","C","A","B",1,1],
        ]
        
        String pedString = pedData*.join('\t').join('\n')+'\n'

        def ped = Pedigrees.parse(new StringReader(pedString))
        
        ped.renameFamily("F","G")
        
        assert ped.families["G"] != null
        assert ped.subjects["B"].id == "G"
    }
}
