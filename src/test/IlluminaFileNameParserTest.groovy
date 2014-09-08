import static org.junit.Assert.*;

import org.junit.Test;


class IlluminaFileNameParserTest {
    
    def parser = new IlluminaFileNameParser()

    @Test
    public void test() {
        
        def info = parser.parse("L10248_S12_L001_R2_001.fastq.gz")
        println info
        
        info = parser.parse("GHS025_D1EH2ACXX_ATCACG_L001_R1.fastq.bam")
        println info
        
        info = parser.parse("008_11_H14NGADXX_L1.1.fastq.trim.atrim.reorder.bam")
        println info
        
        info = parser.parse("../../data/11D14124_M_H14NGADXX_L1.1.fastq.gz")
        assert info.sample == "11D14124_M"
        println info
        
        info = parser.parse("../../data/11D14124_H14NGADXX_L1.1.fastq.gz")
        assert info.sample == "11D14124"
        println info
        
        info = parser.parse("025_11_H14NGADXX_L1_R1.fastq.gz")
        assert info.sample == "025_11"
        println info
         
        
         
    }

}
