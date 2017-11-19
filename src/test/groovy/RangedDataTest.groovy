import static org.junit.Assert.*;

import org.junit.Test;

import gngs.Regions

class RangedDataTest {

    public RangedDataTest() {
    }
    
    @Test
    void testParse() {
        RangedData r = new RangedData(new ByteArrayInputStream(
              """
              chr\tstart\tend\tname\tage
              chr1\t100\t120\tone\tsimon
              chr1\t140\t210\ttwo\tfred
              chr1\t190\t250\tthree\tjane
              chr1\t300\t350\tfour\ttom
              """.stripIndent().trim().bytes
            ).newReader(), 0,1,2
        ).load() 
        assert r[0].name == "one"
        assert r[0].age == "simon"
    }
    
    @Test
    void testPreserveData() {
        RangedData r = new RangedData(new ByteArrayInputStream(
              """
              chr\tstart\tend\tname\tage
              chr1\t100\t120\tone\tsimon
              chr1\t140\t210\ttwo\tfred
              chr1\t190\t250\tthree\tjane
              chr1\t300\t350\tfour\ttom
              """.stripIndent().trim().bytes
            ).newReader(), 0,1,2
        ).load()         
        
        def x = r.grep { it.name == "two" }
        
        assert x[0].name == "two"
        
        def y = x as Regions
        assert y[0].name == "two"
    }
}
