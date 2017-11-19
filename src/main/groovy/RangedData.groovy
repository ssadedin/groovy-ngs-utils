import com.xlson.groovycsv.PropertyMapper;

import gngs.GRange
import gngs.Region
import gngs.Regions

import java.util.zip.GZIPInputStream

import org.codehaus.groovy.runtime.StackTraceUtils;

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
 * <p>
 * RangedData extends the Regions class, so it supports all the usual methods for working with
 * genomic ranges. The only difference is that the ranges involved acquire properties 
 * corresponding to the other columns in the input file.
 * 
 * @author ssadedin@gmail.com
 */
class RangedData extends Regions {
    
    Reader source = null
    
    int chrColumn 
    int startColumn
    int endColumn
    String separator='\t'
    
    int genomeZeroOffset=0
    
    List<String> columns

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
        this(getReader(sourceFile), chrColumn, startColumn, endColumn)
    }
    
    public RangedData(Reader reader, int chrColumn, int startColumn, int endColumn) {
        source = reader
        this.chrColumn = chrColumn
        this.startColumn = startColumn
        this.endColumn = endColumn
    }
    
    RangedData load(Map options=[:], Closure c=null) {
        
        int lineNumber = 0
        
        boolean stripChr = options.stripChr ? true : false
        
        // Some data files (looking at you UCSC) are zero-based instead of 1-based
        if(options.zeroBased)
            genomeZeroOffset=1
        
        // Assume columns on first line
        TSV tsv = new TSV(options, source)
        PropertyMapper currentLine
        try {
            for(PropertyMapper line in tsv) {
                currentLine = line
                if(line.values == null) {
                    // I don't know why sometimes at EOF we get a null values?
                    continue
                }
                    
                Region r = parseRegion(line)
                if(stripChr && r.chr.startsWith('chr')) {
                    r.setChr(r.chr.substring(3)) // must call setChr due to expando
                }
                    
                r.range.extra = r
                if(!columns)
                    columns = line.columns*.key
                    
                line.columns.each { String columnName, int index ->
                    if(index != startColumn && index != endColumn && index != chrColumn) {
                        r.setProperty(columnName, line.values[index])
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
            throw new RuntimeException("Failed to parse line $lineNumber: \n\n" + currentLine.values + "\n\n" + exceptionTrace)
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

    protected Region parseRegion(PropertyMapper line) {
        int startPosition = line.values[startColumn].toInteger()+genomeZeroOffset
        int endPosition = line.values[endColumn].toInteger()+genomeZeroOffset
        Region r = new Region(String.valueOf(line.values[chrColumn]),
                        new GRange(startPosition,endPosition,null))
        return r
    }
    
    static getReader(String fileName) {
       fileName.endsWith(".gz") ? new GZIPInputStream(new FileInputStream(fileName)).newReader() : new File(fileName).newReader()  
    }  
}
