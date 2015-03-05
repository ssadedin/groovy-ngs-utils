import com.xlson.groovycsv.PropertyMapper;


class RangedData extends Regions {
    
    Reader source = null
    
    int chrColumn 
    int startColumn
    int endColumn
    int separator='\t'
    
    int genomeZeroOffset=0

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
        this(new File(sourceFile).newReader(), chrColumn, startColumn, endColumn)
    }
    
    public RangedData(Reader reader, int chrColumn, int startColumn, int endColumn) {
        source = reader
        this.chrColumn = chrColumn
        this.startColumn = startColumn
        this.endColumn = endColumn
    }
    
    RangedData load(Map options=[:], Closure c=null) {
        
        // Some data files (looking at you UCSC) are zero-based instead of 1-based
        if(options.zeroBased)
            genomeZeroOffset=1
        
        // Assume columns on first line
        TSV tsv = new TSV(options, source)
        for(PropertyMapper line in tsv) {
            Region r = new Region(line.values[chrColumn], 
                new GRange(line.values[startColumn].toInteger()+genomeZeroOffset,line.values[endColumn].toInteger()+genomeZeroOffset,null))
            r.range.extra = r
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
        }
        return this
    }
}
