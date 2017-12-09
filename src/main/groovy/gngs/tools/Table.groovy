package gngs.tools

import gngs.Cli
import gngs.Utils
import graxxia.CSV
import graxxia.TSV

class Table {

    static void main(String [] args) {
        
        Cli cli = new Cli(usage: "Table <tsv>", header: "Options:")
        cli.with {
            i 'File to display as table', args:1, required: false
            tsv 'Force TSV format'
            csv 'Force CSV format'
            c 'Columns to show', longOpt: 'columns', args:1, required:false
            x 'Columns to exclude', longOpt: 'exclude', args:1, required:false
            n 'Number of rows to show', args:1, required:false
        }
        
        OptionAccessor opts = cli.parse(args)
        if(!opts) {
            System.exit(1)
        }
        
        def file = opts.i
        if(!opts.i && opts.arguments())
            file = opts.arguments()[0]
            
        def data
        if(!file) {
            data = new TSV(System.in.newReader())
        }
        else
        if(opts.tsv || file.endsWith('tsv')) {
            data = new TSV(file).toListMap()
        }
        else
        if(opts.csv || file.endsWith('csv')) {
            data = new CSV(file).toListMap()
        }
        else
        if(file && new File(file).newReader().readLine().tokenize('\t').size()>3) {
            data = new TSV(file).toListMap()
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
        
        Utils.table(data[0]*.key, data.collect { it*.value })
    }
}
