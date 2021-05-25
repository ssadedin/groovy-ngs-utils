package gngs.tools

import gngs.*
import gngs.plot.Bars
import gngs.plot.Line
import gngs.plot.Plot
import graxxia.IntegerStats
import graxxia.Matrix
import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import groovy.util.logging.Log

import static htsjdk.samtools.CigarOperator.S

/**
 * @author Simon Sadedin
 */
@Log
class DeletionPlot extends ToolBase {
    private MultiCov mcov
    private Region region
    private List<String> bamPaths
    private String sampleOutput
    private String covoStats
    
    /**
     * This object accumulates results from each analysis section that runs and is eventually dumped
     * as output into the file provided as the -json argument
     */
    private ConfigObject jsonResults = new ConfigObject()
    
    int splitReadBinSize = 100
    
    SAM bam

    static void main(args) {
        cli('DeletionPlot -bam <sample of interest> <control bams>', args) {
            region 'Region to plot in chromosomal coordinates (chr:start-end)', args:1, required: true
            covo 'Coverage file output', args:1, required: true
            covoStats 'Coverage stats file output', args:1, required: false
            sample 'The sample of interest', args:1, required: true
            sampleOutput 'Use this value instead of actual sample value when writing JSON', args:1, required: false
            covplot  'Name of file to write coverage plot in', args:1, required: true
            srplot 'Name of the file to write the split read plot in', args:1, required: true
            json 'Combined json output file', args:1, required: false
        }
    }

    @Override
    public void run() {
        
        System.properties['java.awt.headless']='true'
        sampleOutput = opts.sampleOutput ?: opts.sample
        bamPaths = opts.arguments()

        String bamPath = bamPaths.find {
            new SAM(it).samples[0] == opts.sample
        }
        
        if(!bamPath) 
            throw new IllegalArgumentException("The sample specified ($opts.sample) is not present in any of the BAM files provided")

        bam = new SAM(bamPath)

        region = new Region(opts.region)

        createCoverageStats()

        plotCoverage()
        
        plotSplitReads() 
        
        if(opts.json)
            saveJson(opts.json)
    }
    
    void saveJson(String path) {
        log.info "Saving json to $path"
        Utils.writer(path).withWriter { w ->
            w << JsonOutput.prettyPrint(JsonOutput.toJson(jsonResults))
        }
    }

    void createCoverageStats() {
        log.info "Creating coverage stats"

        mcov = new MultiCov()
        Map multicovOpts = [
                rel              : false,
                std              : true,
                '2pass'          : true,
                headers          : true,
                w                : 10,
                subs             : 10,
                o                : opts.covo,
                L                : opts.region,
                covoStats        : opts.covoStats,
                covoSampleOutputs: [("${opts.sample}".toString()): "${opts.sampleOutput}"],
                statsMaxPctValue : 50000,
                arguments        : bamPaths
        ]

        mcov.opts = new CliOptions(overrides:multicovOpts)

        mcov.run()
    }

    void plotCoverage() {
        
        mcov = new MultiCov()
        Map multicovOpts = [
            rel: true,
            std: true,
            '2pass' : true,
            headers: true,
            w : 10,
            subs: 10,
            o : opts.covo,
            L : opts.region,
            arguments: bamPaths
        ]
        
        mcov.opts = new CliOptions(overrides:multicovOpts)

        mcov.run()
        
       
        Matrix covs = Matrix.load(this.opts.covo, r:true)
        
        List covValues = covs[opts.sample] as List
        
        double maxCovValue = covValues.max()
        double plotRangeY = Math.min(1.5, maxCovValue)
        
        List posValues = covs.pos.collect { it }

        Plot covPlot = 
            new Plot(
                title: "Normalised Coverage for ${opts.sample} against ${bamPaths.size()-1} controls" ,
                yBound:[0,plotRangeY], 
                xBound:[region.from, region.to],
                xLabel: 'Position'
             ) << 
            new Line(x: posValues , y: covValues) 

        covPlot.save(this.opts.covplot)
        
        jsonResults[opts.sampleOutput].coverage = [posValues, covValues].transpose().collect { pos, cov ->
            [
                pos: pos, 
                cov: cov
            ]
        }
    }
    
    void plotSplitReads() {
        
        double testSampleMean = this.mcov.means[opts.sample]
        
        Map<String, Matrix> allSplitReads = bamPaths.collect {
            new SAM(it)
        }
        .grep {
            it.samples[0] != opts.sample
        }
        .collectEntries {
            def splitReads = countSplitReads(it)
            double sampleMean = mcov.means[it.samples[0]]
            def downsampled_counts = splitReads.col('count').collate(splitReadBinSize).collect { testSampleMean * it.sum() / sampleMean }
            [it.samples[0], downsampled_counts]
        }
        
        Matrix controls = new Matrix(allSplitReads)
        
        List controlStats = controls.collect {
            new IntegerStats(20000,it)
        }
        
        log.info "Calculated control matrix: \n$controls"
        
        def splitReads = countSplitReads(bam)
       
        final int halfSplitReadBinSize = (int)(splitReadBinSize / 2)

        def downsampled_pos = splitReads.pos.collate(splitReadBinSize).collect { it[halfSplitReadBinSize]}
        
        def downsampled_counts = splitReads.col('count').collate(splitReadBinSize).collect { it.sum() }
        
        def normalised_counts = [controlStats, downsampled_counts].transpose().collect { IntegerStats positionStats, sampleCount ->
            Math.max(0,sampleCount - positionStats.getPercentile(50))
        }
        
        def median_sr = new IntegerStats(10000, normalised_counts).getPercentile(50)
        
        log.info "Calculated median split read count = $median_sr"
        
        def median_subtracted_norm_counts = normalised_counts.collect { Math.max(0, it - median_sr) }
        
        def sr_ds = new Matrix(pos: downsampled_pos, count: median_subtracted_norm_counts)
        
        int maxCount = normalised_counts.max()
        
        Plot plot = new Plot(
            title: "Split Read Counts for ${opts.sample}",
            xBound:[region.from, region.to],
            xLabel: 'Position',
            yBound: [0, Math.max(100,maxCount)]
        ) <<
            new Bars(x:sr_ds.pos as List, y: sr_ds.count as List)
        
        jsonResults[opts.sampleOutput].splitReads = sr_ds.collect {
            [
                pos: pos,
                count: count
            ]
        }.grep {
            it.pos != null
        }

        plot.save(opts.srplot)
    }
    
    final int minSplitReadLength = 7
    
    @CompileStatic
    Matrix countSplitReads(SAM bam) {
        
        log.info "Count split reads in $bam.samFile"
        
        final int regionSize = (int)region.size()

        final int [] srCounts = new int[region.size()+1]
        
        bam.withIterator(region) { i ->
            i.each { r ->
       
                def c0 = r.cigar.cigarElements[0]
                if(!c0)
                    return
                def ce = r.cigar.cigarElements[-1]
        
                if(c0.operator == S && c0.length > minSplitReadLength) {
                    srCounts[r.alignmentStart-region.from]++
                }
                if(ce.operator == S && ce.length > minSplitReadLength) {
                    int srPos = r.alignmentEnd
                    int index = srPos-region.from
                    if(index<regionSize)
                        srCounts[index]++
                }
            }
        }
        
        Matrix splitReads = new Matrix((Map)[pos:region.from..region.to, count: srCounts[0..regionSize]])
        
        return splitReads
        
    }
}
