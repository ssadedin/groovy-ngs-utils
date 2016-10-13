import groovy.lang.Closure;
import groovy.time.TimeCategory;
import groovy.transform.CompileStatic;
import java.text.NumberFormat
import java.util.logging.*
import java.text.*

/**
 * The default Java log former uses a format that is too verbose, so
 * replace it with something more compact.
 */
public class SimpleLogFormatter extends Formatter {
    
    private static final String lineSep = System.getProperty("line.separator");
    
    /**
     * A Custom format implementation that is designed for brevity.
     */
    public String format(LogRecord record) {
        
        DateFormat format = new SimpleDateFormat("h:mm:ss");
    
        String loggerName = record.getLoggerName();
        if(loggerName == null) {
            loggerName = "root";
        }
        StringBuilder output = new StringBuilder()
            .append(loggerName)
            .append("\t[")
            .append(record.threadID)
            .append("]\t")
            .append(record.getLevel()).append("\t|")
            .append(format.format(new Date(record.getMillis())))
            .append(' ')
            .append(record.getMessage()).append(' ')
            .append(lineSep);
            
        if(record.getThrown()!=null) {
            StringWriter w = new StringWriter()
            record.getThrown().printStackTrace(new PrintWriter(w))
            output.append("Exception:\n" + w.toString())
        }    
            
        return output.toString();
    }
}



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

	@CompileStatic
	static int[] array(int... values) {
		return values
	}
    
    @CompileStatic
	static int max(int[] values) {
		int max = Integer.MIN_VALUE;
		for(int value : values) {
				if(value > max)
						max = value;
		}
		return max;
    }
    
    @CompileStatic
	static int min(int[] values) {
		int min = Integer.MAX_VALUE;
		for(int value : values) {
				if(value < min)
						min = value;
		}
		return min;
    }
    
    final static NumberFormat humanNumberFormat = NumberFormat.numberInstance
    
    static {
        humanNumberFormat.maximumFractionDigits=1
        humanNumberFormat.minimumFractionDigits=0
    }
    
    static String humanSize(Number number, units=['bp','kb','Mb','Gb']) {
        double value = number.toDouble()
        for(unit in units) {
            if(value < 1000)
                return humanNumberFormat.format(value) + unit
            value = value / 1000
        }
    }
    
    /**
     * Set up Java logging with a nice, simple, reasonable format
     */
    public static void configureSimpleLogging(Level level = Level.INFO) {
        ConsoleHandler console = new ConsoleHandler()
        console.setFormatter(new SimpleLogFormatter())
        console.setLevel(level)
        Logger log = Logger.getLogger("dummy")
        def parentLog = log.getParent()
        parentLog.getHandlers().each { parentLog.removeHandler(it) }
        log.getParent().addHandler(console)
    }
}
