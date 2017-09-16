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
    
    static timeMs(Closure c) {
        long startTimeMs = System.currentTimeMillis()
        try {
            c()
        }
        finally {
            return System.currentTimeMillis() - startTimeMs
        }
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
        
        // Disable logging from groovy sql
        Logger.getLogger("groovy.sql.Sql").useParentHandlers = false
    }
    
    static String table(Map options = [:], List<String> headers, List<List> rows) {
        
        String indent = options.indent ? (" " * options.indent) : ""
        
        // Create formatters
        Map formatters = options.get('format',[:])
        headers.each { h ->
            if(!formatters[h])
                formatters[h] = { String.valueOf(it) }
            else 
            if(formatters[h] instanceof Closure) {
                // just let it be - it will be called and expected to return the value
            }
            else { // Assume it is a String.format specifier (sprintf style)
                String spec = formatters[h]
                if(spec == "timespan") {
                    formatters[h] = { times ->
                        TimeCategory.minus(times[1],times[0]).toString().replaceAll(TRIM_SECONDS, '').replaceAll(TRIM_ZEROS,' seconds')
                    }
                }
                else {
                    formatters[h] = { val ->
                        String.format(spec, val)
                    }
                }
            }
        }
        
        // Create renderers
        Map renderers = options.get('render',[:])
        headers.each { hd ->
            if(!renderers[hd]) {
                renderers[hd]  = { val, width  -> print val.padRight(width) }
            }
        }
        
        // Find the width of each column
        Map<String,Integer> columnWidths = [:]
        if(rows) {
            headers.eachWithIndex { hd, i ->
                Object widestRow = rows.max { row -> formatters[hd](row[i]).size() }
                columnWidths[hd] = Math.max(hd.size(), formatters[hd](widestRow[i]).size())
            }
        }
        else {
            headers.each { columnWidths[it] = it.size() }
        }
            
        // Now render the table
        String header = headers.collect { hd -> hd.center(columnWidths[hd]) }.join(" | ")
        println indent + header
        println indent + ("-" * header.size())
        
        rows.each { row ->
            int i=0
            headers.each { hd -> 
                if(i!=0)
                    print(" | ");
                else
                    print(indent)
                    
                 renderers[hd](formatters[hd](row[i++]), columnWidths[hd])
            }
            println ""
        }
    }
    
    public static int getNumberOfOpenFiles() {
        java.lang.management.OperatingSystemMXBean os = java.lang.management.ManagementFactory.getOperatingSystemMXBean();
        if(os instanceof com.sun.management.UnixOperatingSystemMXBean){
            return ((com.sun.management.UnixOperatingSystemMXBean) os).getOpenFileDescriptorCount();
        }
    }
}
