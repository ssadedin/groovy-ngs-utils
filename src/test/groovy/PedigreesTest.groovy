import static org.junit.Assert.*;

import org.junit.Test;


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
}
