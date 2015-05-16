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
public class CoverageStats extends graxxia.IntegerStats {
    
    private static final long serialVersionUID = 1L;

    public CoverageStats(int maxPercentileValue, InputStream inStream)
            throws IOException {
        super(maxPercentileValue, inStream);
    }

    public CoverageStats(int maxPercentileValue, Iterable covValues) {
        super(maxPercentileValue, covValues);
    }

    public CoverageStats(int maxPercentileValue) {
        super(maxPercentileValue);
    }
}
