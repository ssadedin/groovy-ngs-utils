package gngs.tools

import com.opencsv.CSVWriter
import gngs.Cli
import gngs.RefGenes
import gngs.Region
import gngs.Utils
import graxxia.CSV
import graxxia.TSV

class Table {
    
    static Map<String,List<String>> AUTO_COLUMN_NAMES = [
        bed: ['Chr', 'Start','End','Id'],
        ped: ['Family','Sample','Father','Mother','Sex','Phenotype'],
        vcf: { path -> ['Chr', 'Pos', 'Ref', 'Alt'] + new gngs.VCF(path).samples }
    ]

    static void main(String [] args) {
        
        Cli cli = new Cli(usage: "Table <tsv>", header: "Options:")
        cli.with {
            i 'File to display as table', args:1, required: false
            t 'Title for table', longOpt: 'title', args:1, required: false
            tsv 'Force TSV format'
            csv 'Force CSV format'
            c 'Columns to show', longOpt: 'columns', args:1, required:false
            x 'Columns to exclude', longOpt: 'exclude', args:1, required:false
            filter 'Row filter: groovy expression to limit rows', args:Cli.UNLIMITED, required:false
            sort 'Sort by given column (specify title)', args: 1, required:false
            gene 'Add a column with an HGNC gene symbol, based on interpreting the first 2 column BED style coordinates', required:false
            n 'Number of rows to show', args:1, required:false
            ofmt 'Output format: csv,tsv,txt default is text', args:1, required: false
            // multi 'If there are empty lines, treat as multiple tables' // todo
            h 'Specify headers. If specified first row of data is treated as data', longOpt: 'headers', args:1, required:false
        }
        
        OptionAccessor opts = cli.parse(args)
        if(!opts) {
            System.exit(1)
        }
        
        List files = []
        if(opts.i) 
            files = [opts.i]
        else
        if(opts.arguments())
            files = opts.arguments() as List
            
        String fileExt = files[0] ? files[0].replaceAll(/^.*\./,'') : null

        Map readOptions = [:]
        if(opts.h) {
            readOptions.columnNames = opts.h.tokenize(',')*.trim()
        }
        else
        if(fileExt in AUTO_COLUMN_NAMES) {
            def autoMapper = AUTO_COLUMN_NAMES[fileExt]
            if(autoMapper instanceof List)
                readOptions.columnNames = autoMapper
            else
            if(autoMapper instanceof Closure)
                readOptions.columnNames = autoMapper.call(files[0])
            else
                assert false
        }
        
        boolean multiFile = files?.size()>1
                        
        def data
        if(!files) {
            data = new TSV(readOptions,System.in.newReader()).toListMap()
        }
        else
        if(opts.tsv || files[0].endsWith('tsv')) {
            data = files.collect { f -> new TSV(readOptions,f).toListMap().collect { (multiFile?[File: f]:[:]) + it }}.sum()
        }
        else
        if(opts.csv || files[0].endsWith('csv')) {
            data = files.collect { f -> new CSV(readOptions,f).toListMap().collect { (multiFile?[File: f]:[:]) + it }}.sum()
        }
        else
        if(files[0] && new File(files[0]).newReader().readLine().tokenize('\t').size()>3) {
            data = files.collect { f -> new TSV(readOptions,f).toListMap().collect { (multiFile?[File: f]:[:]) + it }}.sum()
        }
        else {
            System.err.println()
            System.err.println "Unable to determine format!"
            System.err.println()
            System.exit(0)
        }
        
        if(opts.n) {
            data = data.take(opts.n.toInteger())
        }
       
        if(opts.filters) {
            data = data.grep { row ->
                opts.filters.every { f ->
                    Eval.x(row, f)
                }
            }
        }
        
        if(opts.sort) {
            String sortField = opts.sort
            Integer multiplier = null
            if(sortField.startsWith('-')) {
                multiplier = -1
                sortField = sortField.substring(1)
            }
            data = data.sort { row ->
                def x = row[sortField]
                if(multiplier != null) {
                    x = x * multiplier
                }
                return x
            }
        }
        
        if(opts.c) {
            Set<String> columns = opts.c.tokenize(',')*.trim() as Set
            data = data.collect { row ->
                row.grep { it.key in columns }.collectEntries()
            }
        }
        
        if(opts.x) {
            Set<String> columns = opts.x.tokenize(',')*.trim() as Set
            data = data.collect { row ->
                row.grep { !(it.key in columns) }.collectEntries()
            }            
        }
        
        String outputFormat = 'txt'
        if(opts.ofmt)
            outputFormat = opts.ofmt
            
        if(opts.gene) {
            annotateGenes(data)
        }
        
        def tableOptions = [:]
        if(opts.title) {
            tableOptions.title = opts.title
        }
        
        if(outputFormat == 'txt') {
            Utils.table(tableOptions, data[0]*.key, data.collect { it*.value })
        }
        else
        if(outputFormat == 'csv') {
            writeCSV(data)
        }
        else
        if(outputFormat == 'tsv') {
            writeTSV(data)
        }
    }
    
    static void writeTSV(List<Map> data) {
        System.out.withWriter { Writer w ->
            w.println(data[0]*.key.join('\t'))
            for(Map row in data) {
                w.println(row*.value.join('\t'))
            }
        }    
    }
    
    static void writeCSV(List<Map> data) {
        writeData(data,',')
    }
    
    static void writeData(List<Map> data, String sep) {
        System.out.withWriter { Writer w ->
            CSVWriter csv = new CSVWriter(w, sep as char)
            csv.writeNext(data[0]*.key as String[])
            for(Map row in data) {
                csv.writeNext(row*.value.collect { String.valueOf(it) } as String [])
            }
        }
    }
    
    static void annotateGenes(List<Map> data) {
        RefGenes refgenes = RefGenes.download()
        data.each { Map row ->
            List values = row*.value
            row.Gene = refgenes.getGenes(new Region(values[0]  , values[1].toInteger(), values[2].toInteger())).join(',')
        }
    }
}
