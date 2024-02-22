package gngs.tools

import gngs.ToolBase
import gngs.plot.Plot
import graxxia.IntegerStats
import groovy.transform.CompileStatic
import groovy.util.logging.Log
import groovyx.gpars.GParsPool
import gngs.*
import graxxia.*

@Log
class InsertSizePlot extends ToolBase {

    @Override
    public void run() {
        
        List<SAM> bams = opts.arguments().collect { path ->
            log.info "Computing insert size statistics for $path" 
            new SAM(path)
        }
        
        int maxInferredInsertSize = opts['maxInferredInsertSize']?:5000
        
        Map<String, IntegerStats> allStats = GParsPool.withPool(opts.t?:1) {
            bams.collectParallel {
                [ it.samples[0], computeInsertSizeStats(it, maxInferredInsertSize) ]
            }.collectEntries()
        }
        
        int step = opts.step ?: 10 
        int max = opts.max ?: 1000
        
        def points = (1..max).step(step)
        
        Plot plot = 
            new Plot(
                title: 'Insert Size Distribution Plot',
                legendLocation: 'east',
                xLabel: 'Insert Size', 
                yLabel: 'Fequency'
            ) 
        
        allStats.each { sample, stats ->
            plot << new gngs.plot.Line(x: points, y: (stats.values as List).collate(step).collect { Stats.mean(it) } , displayName: sample)
        }
        
        
        plot.save(opts.o)
        
        log.info "Saved $opts.o"
        
        def results = allStats.collect { sample, stats ->
                [
                    Sample: sample, 
                    Mean: stats.mean,
                    Median: stats.median,
                    StdDev: stats.standardDeviation
                ]
            }

        if(opts.otsv) {
            graxxia.TSV.save(results, opts.otsv)
        }
        Utils.table(results)
    }
    
    @CompileStatic
    IntegerStats computeInsertSizeStats(SAM bam, int maxInferredInsertSize) {
        IntegerStats s = new IntegerStats(maxInferredInsertSize)
        
        log.info "Process $bam.samples"
        
        if(opts['chr']) {
            bam.withIterator(new Region("${opts['chr']}:0-0")) {
                it.each {
                    if(it.properPairFlag && it.firstOfPairFlag && it.inferredInsertSize < maxInferredInsertSize) s.addValue(Math.abs(it.inferredInsertSize))
                }
            }
        }
        else {
            bam.withIterator {
                it.each {
                    if(it.properPairFlag && it.firstOfPairFlag && it.inferredInsertSize < maxInferredInsertSize) s.addValue(Math.abs(it.inferredInsertSize))
                }
            }
        }
        return s
    }
    
    static void main(String [] args) {
        cli('<bam> [<bam>]...', args) {
            t 'Number of threads to use', args:1, type: Integer
            chr 'Limit to chromsome <chr>', args:1
            step 'Granularity of plot (distance between points plotted)', args:1, type: Integer
            max 'Max insert size to show on plot', args:1, type: Integer
            otsv 'Output high level statistics in tsv form', args:1, type: String
            o 'Output file name', args:1, type: String
        }
    }

}
