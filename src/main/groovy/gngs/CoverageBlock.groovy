package gngs

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics

import groovy.lang.IntRange
import groovy.transform.CompileStatic
import groovy.transform.ToString

/**
 * A run of contiguous bases within a band of coverage depth
 * (eg: a run of low coverage bases).
 * <p>
 * Note that the start and end values are <b>inclusive</b>. That is,
 * each base including both start and end values is below coverage
 * threshold.
 *
 * @author Simon Sadedin
 */
@CompileStatic
@ToString(includeNames=true,excludes=['stats'])
class CoverageBlock implements IRegion {
    String region
    String chr
    int start
    int end
    String id
    List<Map> annotations
    
    DescriptiveStatistics stats = new DescriptiveStatistics()
    
    Object asType(Class clazz) {
        if(clazz == Region) {
            return new Region(chr, new GRange(start,end, this))
        }
    }

    @Override
    public IntRange getRange() {
        return start..end
    }
    
    int size() {
        return end - start + 1
    }
}

