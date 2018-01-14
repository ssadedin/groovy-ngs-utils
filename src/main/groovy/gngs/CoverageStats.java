/*
 *  Groovy NGS Utils - Some simple utilites for processing Next Generation Sequencing data.
 *
 *  Copyright (C) 2018 Simon Sadedin, ssadedin<at>gmail.com
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
package gngs;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.StringTokenizer;
import java.util.zip.GZIPInputStream;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import groovy.transform.CompileStatic;

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
    
    private double totalBases = -1d;

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
    
    public double fractionAbove(int coverageLevel) {
        if(totalBases<0d) {
            long total = 0L;
            for(int i=0; i<this.values.length; ++i) {
                total += values[i];
            }
            totalBases = (double)total;
        }
        
        long total = 0;
        for(int i=coverageLevel+1; i<values.length; ++i) {
            total += values[i];
        }
        return total / totalBases;
//        this.values[x..-1].sum() / totalBases 
    }
    
    @CompileStatic
    static CoverageStats from(InputStream stream) throws IOException {
        CoverageStats result = new CoverageStats(1000);
        BufferedReader reader = new BufferedReader(new InputStreamReader(stream));
        try {
            String line = reader.readLine();
            while(line != null) {
                String covValue = line.substring(line.lastIndexOf('\t')+1);
                int intValue = Integer.parseInt(covValue);
                result.addValue(intValue);
                line = reader.readLine();
            }
        }
        finally {
            try {
                reader.close();
            }
            catch(Exception e) {
                // nothing
            }
        }
        return result;
    }        
}
