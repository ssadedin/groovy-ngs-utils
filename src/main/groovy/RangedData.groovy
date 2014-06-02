import com.xlson.groovycsv.PropertyMapper;


class RangedData extends Regions {
    
    Reader source = null
    
    int chrColumn 
    int startColumn
    int endColumn
    int separator='\t'

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
        // Assume columns on first line
        TSV tsv = new TSV(options, source)
        for(PropertyMapper line in tsv) {
            Region r = new Region(line.values[chrColumn], new GRange(line.values[startColumn.toInteger()],line.values[endColumn.toInteger()],null))
            r.range.extra = r
            line.columns.each { String columnName, int index ->
                if(index != startColumn && index != endColumn && index != chrColumn) {
//                    println "Setting $columnName => " + line.values[index]
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
