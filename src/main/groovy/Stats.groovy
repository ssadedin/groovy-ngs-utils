import groovy.transform.CompileStatic;

import java.text.ParseException;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

class Stats extends DescriptiveStatistics {

    public Stats() {
    }
    
    void leftShift(value) {
        this.addValue(value)
    }
    
    static Stats from(Collection values, Closure c=null) {
        Stats s = new Stats()
        values.each {
              s.addValue(c==null? it : c(it))
        }
        return s
    }
    
    @CompileStatic
    static mean() {
        SummaryStatistics s = new SummaryStatistics()
        int errorCount = 0;
        System.in.eachLine { String line ->
            try {
                s.addValue(Double.parseDouble(line.trim()))
            }
            catch(ParseException e) {
                ++errorCount
            }
        }
        if(errorCount>0)
            System.err.println "WARNING: $errorCount lines could not be parsed as numbers"
            
        return s.mean
    }
    
    static Double mean(Iterable iterable) {
        return summary(iterable.iterator()).mean
    }
    
    static Double mean(Iterator i) {
        return summary(i).mean
    }
     
    static SummaryStatistics summary(Iterable iterable) {
        summary(iterable.iterator())
    }
    
    static summary(Iterator i) {
        summary { i.next() }
    }
    
    static Double mean(Closure c) {
        summary(c).mean
    }
    
    @CompileStatic
    static SummaryStatistics summary(double [] values) {
        int i=0;
        summary { values[i++] }
    }
    
    /**
     * Convenience method for returning the median of values read from stdin
     * using a default maximum value of 10,000.
     * 
     * @return
     */
    static median() {
        percentile().getPercentile(50)
    }
    
     static median(int max, Closure c) {
        percentile(max,c).getPercentile(50)
    }
    
    static median(Closure c) {
        percentile(10000,c).getPercentile(50)
    }
    
     static percentile(int max) {
        percentile(max, System.in)
    }
    
    static percentile() {
        percentile(10000, System.in)
    }
     
    /**
     * Return a CoverageStats object by reading lines from the
     * given input stream until no more lines are left.
     * 
     * @param max   See #percentile(max,Closure)
     * @param i     input stream
     * @return  a CoverageStats object
     */
    @CompileStatic
    static percentile(int max, InputStream i) {
        i.withReader { Reader r ->
            percentile(10000) {
                r.readLine()?.trim()
            }
        }
    }
    
    /**
     * Calculate Percentile object by accepting values from
     * the given closure until it either:
     * 
     * <li>returns null
     * <li>throws NoSuchElementException
     * <li>throws ArrayIndexOutOfBounds exception
     * 
     * @param max   An estimate of the maximum possible value that the percentile
     *              values that are requested will have. If the actual value is 
     *              above this level then an exception will be thrown when the
     *              value is requested
     * @param c
     * @return
     */
    @CompileStatic
    static percentile(int max, Closure c) {
        CoverageStats p = new CoverageStats(max)
        def value = c()
        try {
            while(value != null) {
                if(value instanceof Integer) {
                    p.addValue(value)
                }
                else
                if(value instanceof Double) {
                    p.addValue(((Double)value).toInteger())
                }
                else
                if(value instanceof Float) {
                    p.addValue(((Float)value).toInteger())
                }
                else
                if(value instanceof String) {
                    p.addValue(((String)value).toInteger())
                }
                else {
                    p.addValue(String.valueOf(value).toInteger())
                }
                value = c()
            }
        }
        catch(NoSuchElementException e) {
            // Do nothing!
        }
        catch(ArrayIndexOutOfBoundsException e) {
            // Do nothing!
        }
        return p
    }
    
    /**
     * Return a SummaryStatistics object obtained by executing the
     * given closure repeatedly until it either 
     * 
     * <li>returns null
     * <li>throws NoSuchElementException
     * <li>throws ArrayIndexOutOfBounds exception
     * 
     * @param c
     * @return {@link SummaryStatistics} object
     */
    @CompileStatic
    static SummaryStatistics summary(Closure c) {
        SummaryStatistics s = new SummaryStatistics()
        def value = c()
        try {
            while(value != null) {
                if(value instanceof Double) {
                    s.addValue((Double)value)
                }
                else
                if(value instanceof Float) {
                    s.addValue(((Float)value).toDouble())
                }            
                else
                if(value instanceof Integer) {
                    s.addValue(((Integer)value).toDouble())
                }                    
                else
                if(value instanceof String) {
                    s.addValue(((String)value).toDouble())
                }
                else {
                    s.addValue(String.valueOf(value).toDouble())
                }
                value = c()
            }
        }
        catch(NoSuchElementException e) {
            // Do nothing!
        }
        catch(ArrayIndexOutOfBoundsException e) {
            // Do nothing!
        }
        return s
    }
}
