

import org.junit.Test;

class GATKIntervalSummaryTest {
    
    @Test
    void testReadIntervalSummary() {
        
        GATKIntervalSummary summary = 
            new GATKIntervalSummary("src/test/data/test.sample_interval_summary").load()
        
        println "Loaded ${summary.numberOfRanges} intervals"
        
        
        def cov = summary[0].average_coverage
        
        println "Coverage of first interval = " + cov + " (type="+cov.class.name+")"
        
        println "The average coverage depth = " + (summary*.average_coverage.sum() / summary.numberOfRanges)
        
    }

} 
