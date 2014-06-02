import groovy.lang.Closure;
import groovy.time.TimeCategory;


/**
 * Miscellaneous utilities that I couldn't think to put anywhere else
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class Utils {
    static time(String desc, Closure c) {
        System.err.println((" Starting " + desc + " ").center(80, "="))
        Date startTime = new Date()
        Date endTime = startTime
        try {
            return c()
        }
        finally {
            endTime = new Date()
            System.err.println(("$desc executed in " + TimeCategory.minus(endTime,startTime)).center(80,"="))
        }
        // return endTime.time - startTime.time
    }
}
