import static org.junit.Assert.*;

import org.junit.Test;


class SampleInfoTest {
	
//	static String TEST_FILE_NAME = "test_sample_info.tsv"
	static String TEST_FILE_NAME = "samples_AGRF_002.txt"

	@Test
	public void testParse() {
        
		if(!new File("tests/$TEST_FILE_NAME").exists()) {
			println "Skipping test parse because file does not exist. Please download and place '$TEST_FILE_NAME' in the tests local directory."
			return
		}
		
		Map samples = SampleInfo.parse_sample_info("tests/$TEST_FILE_NAME")
		
		SampleInfo s = samples["010100101"]
		
		println "Checking sample $s"
		
		assert s != null
		
		assert s.batch == "002"
		assert s.target == "CS"
		assert s.files.fastq.find { it.startsWith("../data/Sample1_Lab2_Batch1_")}
		assert s.sex == Sex.MALE
		assert s.ethnicity == Ethnicity.EUROPEAN
		assert s.consanguinity == Consanguinity.SUSPECTED
		assert s.sampleType == SampleType.NORMAL
		assert s.dnaConcentrationNg == 30
		assert s.dnaQuantity == 550
		assert s.dnaQuality == 1.8f
		assert s.dnaQuality == 1.8f
		assert s.sequencingDates[0].toCalendar().get(Calendar.DAY_OF_MONTH) == 10
		assert s.sequencingContact == "sequencing@melbournegenomics.org.au"
		assert s.analysisContact == "pipeline@melbournegenomics.org.au"
		assert s.machineIds[1] == "700890"
		
		s.validate()
	}
	
	@Test
	void testNoOptionalFields() {
		
		Map samples = SampleInfo.parse_sample_info("tests/samples.simple.txt")
		
		SampleInfo s = samples["000000000"]
		
		println "Parsed sample $s"
		assert s.batch == "000"
		
	}
	
	@Test
	void testParseDate() {
		SampleInfo s = new SampleInfo()
		Calendar d = s.parseDate("20140320").toCalendar()
		assert d.get(Calendar.DAY_OF_MONTH) == 20
		assert d.get(Calendar.YEAR) == 2014
		assert d.get(Calendar.MONTH) == 2
	}
    
    @Test
    void testSampleInfoFromFiles() {
        def samples = SampleInfo.fromFiles(null)
    }
    
    @Test
    void testMgFormat() {
        def sinfo = SampleInfo.parse_mg_sample_info("/Users/simon/work/mg/batches/001/samples.txt")
    }

}
