import static org.junit.Assert.*;

import org.junit.Test;


class TSVTest {

    @Test
    public void test() {
        String testTsv = 
            [ 
               ["foo","cat","5", "10"],
               ["bar","dog","9", "4.2"]
            ]*.join("\t").join("\n").trim().stripIndent()
            
        println "Parsing: $testTsv"
        for(line in new TSV(new StringReader(testTsv), columnNames: ["name","species", "age","weight"])) {
            println "Line = $line"
            println "$line.name is a $line.species and is $line.age years old, weighing $line.weight"
            
            println line.age.class.name
            
            assert line.age instanceof Integer
            assert line.weight instanceof Double
        }
    }
    
    @Test
    public void testGZipInput() {
        new TSV("tests/test.gz").each {
            println "Line: " + it
        }
    }

}
