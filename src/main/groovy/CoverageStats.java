import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Arrays;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

/**
 * An efficient method to calculate percentiles of coverage values
 * that doesn't require holding them all in memory or sorting them.
 * <p>
 * It relies on some limitations regarding coverage values:
 * <li>They are integers
 * <li>We usually know the approximate range, and don't care about it
 *     above a certain value. For example, we are unlikely to observe
 *     values above 1000 and those that do are unlikely to affect the
 *     50'th percentile.
 * 
 * @author simon.sadedin@mcri.edu.au
 */
public class CoverageStats extends SummaryStatistics {
    
    private static final long serialVersionUID = 1L;

    int [] values = null;
    
    int total = 0;
    
    /**
     * 
     * @param maxPercentileValue
     */
    public CoverageStats(int maxPercentileValue) {
        values = new int[maxPercentileValue];
        Arrays.fill(values, 0);
    }
    
    public CoverageStats(int maxPercentileValue, InputStream inStream) throws IOException {
        values = new int[maxPercentileValue];
        BufferedReader r = new BufferedReader(new InputStreamReader(inStream));
        String line = null;
        while((line = r.readLine()) != null) {
            leftShift(line);
        }
    }
     
    void leftShift(Object obj) {
        if(obj instanceof Integer) {
            addValue((Integer)obj);
        }
        else 
        if(obj instanceof String) {
            addValue(Integer.parseInt(String.valueOf(obj).trim()));
        }
        else 
        if(obj instanceof Number) {
            addValue(((Number)obj).intValue());
        }    
     }
    
    /**
     * Count the given coverage value in calculating the median
     * @param coverage
     */
    void addValue(int coverage) {
        if(coverage>=values.length)
            ++values[values.length-1];
        else
            ++values[coverage];
        
        super.addValue(coverage);
        ++total;
    }
    
    /**
     * Return the specified percentile from the observed coverage counts
     * eg: for median, getPercentile(50).
     * 
     * @return the specified percentile, if it is smaller than the max value passed in the
     *         constructor.
     */
    int getPercentile(int percentile) {
        int observationsPassed = 0;
        int lowerValue = -1;
        final int medianIndex = total / (100/percentile);
        for(int i=0; i<values.length; ++i) {
            observationsPassed += values[i];
            if(observationsPassed >= medianIndex) {
                if(total%2 == 0) {
                    // Find the next value and return average of this one and that
                    lowerValue = i;
                    for(int k=i+1; k<values.length; ++k) {
                        if(values[k]>0) {
                            return (lowerValue + k) / 2;
                        }
                    }
                }
                else
                    return i;
            }
        }
        return -1;
    }
    
    public int getMedian() {
        return getPercentile(50);
    }
    
    public int getAt(int percentile) {
        return getPercentile(percentile);
    }
    
    public String toString() {
        return super.toString() + "Median: " + getMedian() + "\n";
    }
}