import groovy.lang.Closure;
import groovy.time.TimeCategory;


/**
 * Miscellaneous utilities that I couldn't think to put anywhere else
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class Utils {
    static time(String desc, Closure c) {
        Date startTime = new Date()
        try {
            c()
        }
        finally {
            Date endTime = new Date()
            System.err.println "$desc executed in " + TimeCategory.minus(endTime,startTime)
        }
    }
}
