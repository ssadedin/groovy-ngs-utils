package gngs
import java.util.zip.GZIPInputStream

import org.codehaus.groovy.runtime.StackTraceUtils;

import com.xlson.groovycsv.PropertyMapper;

import graxxia.DirectReaderFactory
import graxxia.ReaderFactory
import graxxia.StringReaderFactory
import graxxia.TSV
import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.transform.stc.ClosureParams
import groovy.transform.stc.SimpleType

/**
 * RangedData represents a set of genomic regions with data attached. The data is
 * parsed from a tab separated file, 3 columns of which are expected to be the genomic
 * coordinates. The other columns are parsed and accessible as properties directly on
 * the contained Region objects that are loaded.
 * <p>
 * To use RangedData, first create a RangedData object, and then call the load() method:
 * <pre>
 * // my_file.tsv has columns in first line - chr, start, end, depth
 * r = new RangeData("my_file.tsv").load()
 * </pre>
 * This form assumes that the first 3 columns of the file specify a genomic position in BED
 * style representation (chromosome, start, end). These first three columns are extracted and the 
 * remaining columns are sniffed to infer their types, and then added as expandos to the created 
 * Region objects. The result is that you can access them directly as properties.
 * <pre>
 * // Find the 'depth' property of all ranges overlapping chr1:100000-200000
 * println r.grep { it.overlaps("chr1",100000,200000) }.depth
 * </pre>
 * If the file doesn't have column names as the first line then you should specify them yourself:
 * <pre>
 * r = new RangeData("my_file.tsv").load(columnNames:['chr','start','end','depth'])
 * </pre>
 * Note: the 'chr','start' and 'end' columns here arbitrary - they don't affect what is 
 * parsed into the ranges.
 * <p>If you are loading a CSV file, pass the separator to the "load" function:</p>
 * <pre>
 *  r = new RangeData("my_file.csv").load(separator:',')
 * </pre>
 * <p>
 * RangedData extends the Regions class, so it supports all the usual methods for working with
 * genomic ranges. The only difference is that the ranges involved acquire properties 
 * corresponding to the other columns in the input file.
 * 
 * @author ssadedin@gmail.com
 */
class RangedData extends Regions {
    
    ReaderFactory source = null
    
    /**
     * Index of the column containing the reference sequence (or "chromosome")
     */
    int chrColumn 
    
    /**
     * Index of the column containing the start index
     */
    int startColumn
    
    /**
     * Index of the column containing the end index
     */
    int endColumn
    
    /**
     * The separator used between values in the file
     */
    String separator='\t'
    
    /**
     * The starting index for the first base in the genome. Some formats use 1 as the first base,
     * but mostly it is zero.
     */
    int genomeZeroOffset=0
    
    List<String> columns
    
    Closure regionParser = null

    public RangedData() {
    }
    
    /**
     * Default to first 3 columns of file being the genomic range information
     * in form of chr,start,end
     * 
     * @param sourceFile
     */
    public RangedData(String sourceFile) {
        this(sourceFile,0,1,2)
    }
    
    public RangedData(String sourceFile, int chrColumn, int startColumn, int endColumn) {
        this(new StringReaderFactory(source:sourceFile), chrColumn, startColumn, endColumn)
    }
    
    public RangedData(Reader reader, int chrColumn, int startColumn, int endColumn) {
        this(new DirectReaderFactory(reader:reader), chrColumn, startColumn, endColumn)
    }
    
    public RangedData(Reader reader, @ClosureParams(value=SimpleType, options=['com.xlson.groovycsv.PropertyMapper']) Closure regionParser) {
        this(new DirectReaderFactory(reader:reader), -1, -1, -1)
        this.regionParser = regionParser
    }

    public RangedData(ReaderFactory reader, int chrColumn, int startColumn, int endColumn) {
        source = reader
        this.chrColumn = chrColumn
        this.startColumn = startColumn
        this.endColumn = endColumn
    }
    
    @CompileStatic
    public static  Regions loadJSON(def source) {
        List<Map> data = (List<Map>)Utils.reader(source) { Reader r ->
            new JsonSlurper().parse(r)
        }
        
        Regions regions = new Regions()
        for(Map<String,Object> d in data) {
            if(!d.containsKey('chr') || !d.containsKey('start') || !d.containsKey('end'))   
                throw new IllegalArgumentException("No 'chr' column in data : please supply JSON format with chr,start,end attributes")

            Region r = new Region((String)d.chr, (int)d.start, (int)d.end)
            d.each { String k, Object v ->
                if(k != 'chr') {
                    r[k] = v
                }
            }
            
            regions.addRegion(r)
        }
        
        return regions
    }
    
    @CompileStatic
    RangedData load(Map options=[:], Closure c=null) {
        
        int lineNumber = 0
        
        boolean stripChr = options.stripChr ? true : false
        
        // Some data files (looking at you UCSC) are zero-based instead of 1-based
        if(options.zeroBased)
            genomeZeroOffset=1
            
        if(options.regionParser) {
            this.regionParser = (Closure) options.regionParser
        }
        
        // Assume columns on first line
        TSV tsv = new TSV(options, source.newReader())
        PropertyMapper currentLine
        try {
            for(PropertyMapper line in tsv) {
                currentLine = line
                if(line.values == null) {
                    // I don't know why sometimes at EOF we get a null values?
                    continue
                }
                    
                Region r = (regionParser.is(null) ? parseRegion(line) : (Region)regionParser(line))

                if(stripChr && r.chr.startsWith('chr')) {
                    r.setChr(r.chr.substring(3)) // must call setChr due to expando
                }
                    
                ((GRange)r.range).extra = r
                
                final Map<String,Object> cols = (Map<String,Object>)line.columns
                if(!columns) {
                    columns = cols*.key
                }
                
                for(Map.Entry<String, Integer> e : cols.entrySet()) {
                    final String columnName = e.getKey();
                    final int index = e.getValue();
                    if(index != startColumn && index != endColumn && index != chrColumn) {
                        r.setProperty(columnName, ((List)line.values)[index]);
                    }
                }

                if(c != null) {
                    if(c(r)==false)
                        continue
                }
                addRegion(r)
                ++lineNumber
            }
            return this
        }
        catch(Exception e) {
            def exceptionTrace = new StringWriter()
            StackTraceUtils.sanitize(e).printStackTrace(new PrintWriter(exceptionTrace))
            throw new RuntimeException("Failed to parse line $lineNumber: \n\n" + currentLine?.values + "\n\n" + exceptionTrace)
        }
    }
    
    /**
     * Convert to a list of map objects
     * 
     * @return
     */
    List<Map<String,Object>> toListMap() {
        this.collect {  [chr: it.chr, start: it.from, end: it.to] + it.properties }
    }

    @CompileStatic
    protected Region parseRegion(final PropertyMapper line) {
        final List<String> lineValues = (List<String>)line.values
        final int startPosition = Integer.parseInt(lineValues[startColumn])+genomeZeroOffset
        final int endPosition = Integer.parseInt(lineValues[endColumn])+genomeZeroOffset
        
        // Manually construct this because in profiling it is a hotspot where groovy
        // spends time choosing the right constructor dynamically
        Region r = new Region()
        r.chr = lineValues[chrColumn]
        GRange range = new GRange(startPosition, endPosition, null)
        range.extra = r
        r.range = range
        
        return r
    }
    
    static getReader(String fileName) {
       fileName.endsWith(".gz") ? new GZIPInputStream(new FileInputStream(fileName)).newReader() : new File(fileName).newReader()  
    }  
}
