package gngs
import groovy.lang.Closure;
import groovy.time.TimeCategory;
import groovy.transform.CompileStatic;
import groovy.transform.NamedDelegate
import groovy.transform.NamedVariant
import groovy.transform.stc.ClosureParams
import groovy.transform.stc.SimpleType
import htsjdk.samtools.util.BlockCompressedInputStream
import htsjdk.samtools.util.BlockCompressedOutputStream

import java.text.NumberFormat
import java.util.logging.*
import java.util.zip.GZIPInputStream
import java.util.zip.GZIPOutputStream

import org.yaml.snakeyaml.Yaml

import java.io.Writer
import java.nio.file.Files
import java.nio.file.Path
import java.text.*

@CompileStatic
class ExecutedProcess {
    Appendable err
    Appendable out
    int exitValue
}

class TableOptions {
    Integer indent
    Map format
    Map render
    Boolean topborder
    Appendable out
    String border
    String title
    Integer precision
    Double color_threshold
    Map<String,Integer> columnWidths
}


/**
 * Miscellaneous utilities that I couldn't think to put anywhere else
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class Utils {
    
    
    /**
     * Call the given closure and print the time it takes to run
     * <p>
     * Options:
     * <li>log - print to the given logger instead of stderr
     * <li>suppressStartMessage - print to the given logger instead of stderr
     * 
     * @param desc
     * @param c
     * @return  value returned by closure C
     */
    @CompileStatic
    static<T> T time(Map options=[:],String desc, Closure<T> c) {
        
        Closure printMsg 
        if(options.log) 
            printMsg = { ((Logger)options['log']).info(it.toString()) } 
        else
            printMsg = { System.err.println(it) }
        
        if(!options.suppressStartMessage)
            printMsg((" Starting " + desc + " ").center(80, "="))
        Date startTime = new Date()
        Date endTime = startTime
        try {
            return c()
        }
        finally {
            endTime = new Date()
            printMsg((" $desc executed in " + TimeCategory.minus(endTime,startTime)).center(80,"="))
        }
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
    
    final static NumberFormat percNumberFormat = NumberFormat.numberInstance
    static {
        percNumberFormat.maximumFractionDigits=2
        percNumberFormat.minimumFractionDigits=0
    }
    
    /**
     * Convenience method to format a fraction in a reasonable way as a percentage
     * <p>
     * Note: includes the '%' symbol.
     */
    static String perc(BigDecimal value) {
            humanNumberFormat.format(100.0d*value)+"%"
    }    
    
    /**
     * Convenience method to format a fraction in a reasonable way as a percentage
     * <p>
     * Note: includes the '%' symbol.
     */
    static String perc(double value) {
            humanNumberFormat.format(100.0d*value)+"%"
    }
    
    /**
     * Convenience method to format a number in a "human readable" manner, where
     * this means rounding and expressing as millions (m), thousands (k), 
     * billions (g) etc.
     */        
    static String human(Number number) {
        humanBp(number,['','k','m','g','t','p'])
    }
    
    
    /**
     * Convenience method to format a number in a "human readable" manner, where
     * this means rounding and expressing as millions (m), thousands (k), 
     * billions (g) etc.
     * <p>
     * Optionally, custom unit labels can be passed as a list in the second
     * parameter.
     */         
    static String humanBp(Number number, units=['bp','kb','Mb','Gb']) {
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
    
    /**
     * Convenience version of {@link #table} to display data in the form of a list of
     * Maps as a table, assuming all the Maps have the same keys which are to be the
     * headers of the table. See {@link TableOptions} for optional parameters.
     * 
     * @param rows  List of map objects representing data, keys of each Map must be identical
     */
    @NamedVariant
    static Writer table(@NamedDelegate TableOptions options=new TableOptions(), List<Map> rows) {
        List headers = rows[0]*.key
        List data = rows.collect { it*.value }
        table(options, headers, data)
    }
    
    final static String BLUE = "\u001B[34m"  // ANSI escape code for blue
    final static String RED = "\u001B[31m"   // ANSI escape code for red
    final static String RESET = "\u001B[0m"  // Reset to default color
    
    /**
     * Pattern to match ansi regex escape sequences with the goal of enabling them to be stripped out
     */
    final static String ANSI_REGEX = ~"\\u001B\\[[0-9;]*m"

    /**
     * A utility to print a table of values in a nice format for 
     * output on a terminal. Columns are aligned, padded, borders
     * drawn etc. The output format is compatible with markdown
     * for downstream re-formatting into documents.
     * <p>
     * By default, the table is printed to stdout. To print it somewhere
     * else, set the <code>out</code> option to a Writer object.
     * <p>
     * Optional parameters: 
     * <li><code>indent</code>: amount to indent table by
     * <li><code>format</code>: Map keyed on column containing custom formatters
     * <li><code>render</code>: Map, keyed on column containing custom renderers
     * <li><code>border</code>: String defining character to use as border
     * <li><code>topborder</code>: if true, a top border will be added
     * <li><code>out</code>: Custom output writer / stream to write results to
     * 
     * @param headers   a list of column names
     * @param rows      a list of lists, where each inner list represents a
     *                  row in the table
     * @return          the writer if one was provided as the <code>out</code> named arg
     */
    @NamedVariant
    static Writer table(@NamedDelegate TableOptions options=new TableOptions(), List<String> headers, List<List> rows) {
        
        String indent = options.indent ? (" " * options.indent) : ""
        if(options.border == 'true') {
            options.border = '|'
        }
        
        def out = options.out ?: System.out
        
        // Create formatters
        Map formatters = options.format?:[:]
        headers.each { h ->
            if(!formatters[h]) {
                formatters[h] = { 
                    String value
                    if(options.precision && it instanceof Number) {
                        value = String.format("%.${options.precision}f", it)
                    }
                    else {
                        value = String.valueOf(it)
                    }
                    
                    if(options.color_threshold && it instanceof Number) {
                       value = (it < options.color_threshold ? BLUE : RED) + value + (options.color_threshold ? RESET : '')
                    }
                    
                    return value
                }
            }
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
        Map renderers = options.render?:[:]
        headers.each { hd ->
            if(!renderers[hd]) {
                renderers[hd]  = { val, width  -> out.print val.padRight(width) }
            }
        }
        
       
        // Find the width of each column
        Map<String,Integer> columnWidths = [:]
        if(rows) {
            headers.eachWithIndex { hd, i ->
                Object widestRow = rows.max { row -> 
                    if(options.color_threshold) {
                        formatters[hd](row[i]).replaceAll(ANSI_REGEX, '').size() 
                    }
                    else
                        formatters[hd](row[i]).size() 

                }
                def widestValue = formatters[hd](widestRow[i])
                columnWidths[hd] = Math.max(hd.size(), options.color_threshold ? widestValue.replaceAll(ANSI_REGEX,'').size() : widestValue.size())
            }
        }
        else {
            headers.each { columnWidths[it] = it.size() }
        }

        if(options.columnWidths) {
            columnWidths = columnWidths + options.columnWidths
        }
        
           
        // Now render the table
        String header = headers.collect { hd -> hd.padRight(columnWidths[hd]) }.join(" | ")
        
        int totalWidth = header.size()
        if(options.title) {
            out.println("-" * (totalWidth))
            if(totalWidth > 160)
                out.println(" " + options.title.padRight((totalWidth-2), " ") + " ")
            else
                out.println(" " + options.title.center((totalWidth-2), " ") + " ")
            out.println("-" * (totalWidth))
        }
 
        String leftBorder = options.border ? options.border + ' ' : ''
        String rightBorder = options.border ? ' ' + options.border : ''

        if(options.topborder) {
            out.println indent + (options.border?' ':'') + ("-" * header.size()) + (options.border?'--':'') 
        }
       
        out.println indent + leftBorder + header + rightBorder
        out.println indent + leftBorder.replaceAll(' ','-') + headers.collect { hd -> '-' *columnWidths[hd] }.join("-|-") + rightBorder.replaceAll(' ','-')
        
        rows.each { row ->
            int i=0
            out.print leftBorder
            headers.each { hd -> 
                if(i!=0)
                    out.print(" | ");
                else
                    out.print(indent)
                    
                renderers[hd](formatters[hd](row[i++]), columnWidths[hd])
            }
            out.println rightBorder
        }
        
        if(options.out) {
            if(out instanceof PrintStream)
                return out.newWriter()
            else
                return out
        }
        else
            return null
    }
    
    public static int getNumberOfOpenFiles() {
        java.lang.management.OperatingSystemMXBean os = java.lang.management.ManagementFactory.getOperatingSystemMXBean();
        if(os instanceof com.sun.management.UnixOperatingSystemMXBean){
            return ((com.sun.management.UnixOperatingSystemMXBean) os).getOpenFileDescriptorCount();
        }
    }
    
    static closeQuietly(obj) {
        if(obj == null)
            return
        try {
            obj.close()
        }
        catch(Exception e) {
            // ignore
        }
    } 
    
    /**
     * Execute the given action repeatedly until it does not throw an exception,
     * or the given number of attempts is exceeded.
     * <p>
     * An optional <code>message</code> argument can be specified to give a 
     * custom message to display with each retry.
     * 
     * @param maxRetries    number of attempts
     * @param action        action to execute
     * @return  the return value from the action
     */
    static Object withRetries(Map options=[:], int maxRetries, Closure action) {
        Logger log = Logger.getLogger('gngs.Utils')
        int count = 0
        long sleepTimeMs = 1000
        if(options.sleepTimeMs != null)
            sleepTimeMs = options.sleepTimeMs.toLong()
            
        while(true) {
             try {
                 return action(count)
             }
             catch(Exception e) {
                 if(count > maxRetries) {
                     log.warning "Exceed max $maxRetries: $options.message"
                     throw e
                 }
                 if(options.message != null)
                     log.info "$options.message: try ${count+1} of $maxRetries failed ($e), retrying ..."
                 else
                     log.info "Try ${count+1} of $maxRetries failed, retrying ..."
             }
                 
             Thread.sleep(sleepTimeMs)
             ++count
             sleepTimeMs = Math.min(10000,sleepTimeMs*2)
         }
     }
     
    /**
     * Creates a Reader from any object that can be reasonably interpreted to represent a file,
     * allowing for the possibility that it might be gzipped or bgzipped.
     * <p>
     * If a {@link groovy.lang.Closure} is provided, then the closure will be called, passing the
     * reader as the first parameter. The return value of this function will be the result of 
     * executing the closure. If, however, a Closure is not provided, the return value will be the
     * created {@link java.io.Reader} object to be used.
     * 
     * @param fileLike
     * @param c Optional {@link groovy.lang.Closure} to call back with the created {@link java.io.Reader}
     * @return  result of evaluated closure if provided, or created {@link java.io.Reader}
     */
    @CompileStatic
    static Object reader(fileLike, @ClosureParams(value=SimpleType, options=['java.io.Reader']) Closure c = null) {
        Reader r = createReader(fileLike)
        if(c != null) {
            Object result
            r.withReader {
                result = c(r)
            }
            return result
        }
        else
            return r
    }
    
    static Reader createReader(fileLike) {
        
        if(fileLike instanceof Reader)
            return fileLike
            
        InputStream stream = Utils.createStream(fileLike)
        
        return stream.newReader()
    }
    
    static InputStream createStream(fileLike) {
       
        boolean gzip = false
        boolean bgzip = false
        
        if(fileLike instanceof String) {
            fileLike = new File(fileLike)
        }
        
        if(fileLike instanceof File) {
            fileLike = fileLike.toPath()
        }
        
        if(fileLike instanceof Path) {
            Path path = fileLike
            if(path.toString().endsWith('.gz'))
                gzip = true
            else
            if(path.toString().endsWith('.bgz'))
                bgzip = true
            fileLike = Files.newInputStream(fileLike)
        }
        
        if(!(fileLike instanceof InputStream))
            throw new IllegalArgumentException("Expected object of type String, File, Path or InputStream, but was passed " + fileLike.class.name)
        
        if(gzip) {
            fileLike = new GZIPInputStream(fileLike, 128*1024)
        }
        else
        if(bgzip) {
            fileLike = new BlockCompressedInputStream(fileLike)
        }
        return fileLike
    } 
    
    /**
     * Executes the given closure, passing as parameters all of the given filenames opened as 
     * Writers using intelligent interpretation of the filenames - eg: if ending in .gz, use gzip, etc.
     * <p>
     * Note: if a file name is null, blank or false, a null will be returned for the corresponding writer.
     * 
     * @param fileNames file names to open
     * @param c         Closure to call
     * @return  result of closure
     */
    static withWriters(List<String> fileNames, Closure c) {
        List<Writer> writers = []
        try {
            for(fileName in fileNames) {
                writers << (fileName ? outputWriter(fileName) : null)
            }
            c(*writers)
        }
        finally {
            for(Writer w : writers) {
                closeQuietly(w)
            }
        }
    }
    
    @CompileStatic
    static Writer writer(File file) {
        return outputWriter(file.path)
    }
  
    @CompileStatic
    static Writer writer(String fileName) {
        return outputWriter(fileName)
    }
    
    @CompileStatic
    static Writer outputWriter(String fileName) {
        new PrintWriter(outputStream(fileName))
    } 
    
    @CompileStatic
    static OutputStream outputStream(String fileName, int bufferSize=1048576) {
        
        if(fileName == '-' || fileName == '/dev/stdout')
            return System.out
        


        if(fileName.endsWith(".bgz"))
            new BufferedOutputStream(new BlockCompressedOutputStream(fileName), bufferSize)
        else
        if(fileName.endsWith(".gz"))
            new GZIPOutputStream(new FileOutputStream(fileName), bufferSize)
        else
            new BufferedOutputStream(new File(fileName).newOutputStream(), bufferSize)
    } 
 
    
    /**
     * Close all the streams associated with the given process
     * ignoring all exceptions.
     * <p>
     * Note: this is necessary because if the streams are not closed
     * this way it seems they can take a while to be closed, even though
     * the process may have ended. If many processes are executed consecutively
     * the file handle limit can be exhausted even though processes are not
     * executed concurrently.
     */
    public static withProcessStreams(Process p, Closure c) {
        try {
            return c()
        }
        finally {
          try { p.inputStream.close() } catch(Throwable t) { }
          try { p.outputStream.close() } catch(Throwable t) { }
          try { p.errorStream.close() } catch(Throwable t) { }
        }
    }

    @CompileStatic
    static ExecutedProcess exec(Map options = [:], String cmd) {
        exec(options, cmd.tokenize() as List, null)
    }
    
    /**
     * Execute the given command and return back a map with the exit code,
     * the standard output, and std err
     * <p>
     * An optional closure will be executed as a delegate of the ProcessBuilder created
     * to allow configuration.
     *
     * @param startCmd  List of objects (will be converted to strings) as args to command
     * @return ExecutedProcess with exitValue, err and out attributes
     */
    @CompileStatic
    static ExecutedProcess exec(Map options = [:], List<Object> startCmd, Closure builder = null) {
        
        List<String> stringified = startCmd*.toString()
        
        ProcessBuilder pb = new ProcessBuilder(stringified)
        if(builder != null) {
            builder.delegate = pb
            builder()
        }
        
        Process p = pb.start()
        ExecutedProcess result = new ExecutedProcess()
        withProcessStreams(p) {
            StringBuilder defaultOut = new StringBuilder()
            Appendable out = (Appendable)options.out ?: defaultOut
            Appendable err = (Appendable)options.err ?: defaultOut
            
            // Note: observed issue with hang here on Broad cluster
            // seems to be related to hang inside OS / NFS call. Maybe use forwarder for this?
            p.waitForProcessOutput(out, err)
            
            result.exitValue = p.waitFor()
            result.err = err
            result.out = out
        }
        
        if((options.throwOnError != false) && (result.exitValue != 0))
            throw new Exception("Command returned exit code ${result.exitValue}: " + stringified.join(" ") + "\n\nOutput: $result.out\n\nStd Err:\n\n$result.err")
            
        return result
    }
    
    @CompileStatic
    static String formatIfNumber(NumberFormat fmt, Object obj) {
        obj instanceof Number ? fmt.format(obj) : String.valueOf(obj)
    }
    
    /**
     * Attempts to parse the given string as a list of numbers and returns the list if it can
     */
    static List toNumberList(String value) {
        List tokens = value.tokenize(',')
        if((tokens.size() > 1) && tokens.every { it.isNumber() }) {
            tokens.collect { convertNum(it) }
        }
        else {
            return null
        }
    }
    
    /**
     * Attempts to convert the given value to a number, returns the number if it can,
     * or the value as a string if it can't.
     */
    static def convertNum(String x) {
        if(x.isInteger())
            return x.toInteger()
        else
        if(x.isDouble())
            return x.toDouble()
        else
            return x
    }
    
    /**
     * Parse a YAML file using default settings
     */
    static Object yaml(def fileLike) {
        Utils.reader(fileLike) {
            return new Yaml().load(it)
        }
    }
}
