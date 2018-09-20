/*
 *  Groovy NGS Utils - Some simple utilites for processing Next Generation Sequencing data.
 *
 *  Copyright (C) 2018 Simon Sadedin, ssadedin<at>gmail.com
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
package gngs.tools

import java.text.NumberFormat

import org.apache.commons.math3.analysis.interpolation.LoessInterpolator
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction

import gngs.*
import gngs.coverage.CoverageCombinerActor
import gngs.coverage.CoveragePrinter
import gngs.coverage.CoverageCalculatorActor
import graxxia.IntegerStats
import graxxia.Matrix
import graxxia.Stats
import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import groovy.util.logging.Log
import groovyx.gpars.GParsPool
import htsjdk.samtools.SAMRecordIterator

/**
 * A command line tool to rapidly calculate per-base coverage across multiple BAM files
 * with extended support for coverage normalisation and variability calculation.
 * <p>
 * MultiCov reads BAM files directly and calculates the coverage over a specified
 * set of regions. It is designed to ensure reading, writing and coverage calculation
 * are all highly parallelised for maximum performance.
 * 
 * @author Simon Sadedin
 */
@Log
class MultiCov extends ToolBase {
    
    Regions scanRegions
    
    List<SAM> bams
    
    NumberFormat fmt = NumberFormat.numberInstance        
    
    
    boolean phase1 = false
    
    static void main(String [] args) {
        cli('-bed <bed file> <bam1> <bam2> ...', args) {
            bed 'BED file describing target regions to analyse coverage for', args:1
            'L' 'A specific region to process, or if a file, opened as a BED file', args:1
            gene 'A specific gene to scan', args:1
            cv 'Output coefficient of variation for each position as column 2', required: false
            rel 'Output values relative to sample mean'
            std 'Standardise output values to mean at each bp'
            stats 'Print statistics about coverage values and variation'
            corr 'Determine corellation of other samples to specified samples (separated by comma)', args:1
            co 'Correlation output file', longOpt: 'correlationOutput', args:Cli.UNLIMITED, required: false
            cvo 'Coefficient of variation output file: write coeffv in tab separated form to this file', args:1, required:false 
            cvj 'Coefficient of variation output file: write coeffv in javascript loadable form to this file', args:1, required:false 
            gcref 'Reference to compute GC content against', args:1, required: false
            '2pass' 'use two phase analysis to get accurate estimates of means'
            o 'Output file to write results to', longOpt: 'output', args:1
            gcprofile 'Write GC profile for each sample to this file in JSON format (requires gcprof)', args:1, required: false 
            covo 'Write statistics about coverage to the given file in js format', args:1
        }
    }
    
    @Override
    public void run() {
        
        long startTimeMs = System.currentTimeMillis()
        
        if(opts.gcprofile) {
            if(!opts.gcref)
                throw new IllegalArgumentException('Please specify the -gcref option to use the -gcprofile option')
        }
        
        fmt.maximumFractionDigits = 2
        
        resolveRegionsToScan()
        
        this.scanRegions = this.scanRegions.reduce().toSorted(new NumericRegionComparator())    
        
        this.bams = opts.arguments().collect { new SAM(it) }
        
        List<String> samples = bams*.samples*.getAt(0)
        log.info "Analysing coverage over ${Utils.humanBp(scanRegions.size())} for ${samples.size()} samples: ${samples.join(',')}"
        
        def output = opts.o ? new File(opts.o) : System.out
        output.withWriter { w ->
            
            
            Map options = [:]
            if(opts.gcref)
                options.gcReference = new FASTA(opts.gcref)
                
            CoveragePrinter printer = new CoveragePrinter(options, w, samples)
            
            if(opts['2pass']) {
                log.info " Executing Phase 1 / 2 to estimate sample means ".center(80,"=")
                CoveragePrinter printer1 = new CoveragePrinter(null, samples)
                this.phase1 = true
                run(printer1)
                
                Map means = samples.collectEntries { sample ->
                    [ sample,  printer1.sampleStats[printer1.samples.indexOf(sample)].mean ]
                }
                printer.setSampleMeans(means)
                
                log.info "Set means to: $means"
                
                log.info "Estimated coverage from ${printer1.sampleStats[0].n} values"
                
                this.phase1 = false
                log.info " Executing Phase 2 / 2 to Compute Coverage Statistics ".center(80,"=")
            }
            
            
            if(opts.cv)
                printer.coeffV = true
            if(opts.rel)
                printer.relative = true
            if(opts.std)
                printer.std = true
            if(opts.corr)
                printer.correlationSamples = opts.corr == '.' ? samples : opts.corr.tokenize(',')
            run(printer)
        }
        
        log.info "Finished in ${Utils.human((System.currentTimeMillis()-startTimeMs)/1000)} seconds"
        
        Thread.sleep(300000)
    } 
    
    /**
     * Combine the various possible options for specifying which region to analyse
     * into a single regions object (see {@link #scanRegions}).
     */
    void resolveRegionsToScan() {
        if(opts.gene) {
            SAM bam = new SAM(opts.arguments()[0])
            String build = bam.sniffGenomeBuild()
            RefGenes geneDb = RefGenes.download(build)
            this.scanRegions = geneDb.getExons(opts.gene)
            log.info "Gene $opts.gene resolved to $scanRegions"
        }
        else
        if(opts.L) {
            if(opts.L.endsWith('.bed') && new File(opts.L).exists())
                this.scanRegions = new BED(opts.L).load()
            else
                this.scanRegions = new Regions([new Region(opts.L)])
        }
        else
        if(opts.bed) {
            log.info "Reading bed file $opts.bed"
            this.scanRegions = new BED(opts.bed).load()
        }
        else {
            System.err.println "Please specify either -L, -bed or -gene to define the regions to process"
            parser.usage()
            System.exit(1)
        }
    }
    
    void run(CoveragePrinter printer) {
        printer.start()
            
        CoverageCombinerActor combiner = new CoverageCombinerActor(consumer: printer, numSamples: bams.size())
        combiner.start()
            
        // By default, estimate mean coverage from the regions to be scanned
        Regions estRegions = scanRegions
        
        // But if a BED file was provided (eg: in addition to -gene) then we can use
        // the whole bed file to estimate coverage
        if(opts.bed)
            estRegions = new BED(opts.bed).load()
        
        GParsPool.withPool(bams.size()) {
            
            // If no means were provided, compute them
            if(printer.sampleMeans == null || printer.sampleMeans.isEmpty()) {
                Map<String, Double> sampleMeans = [:]
                bams.eachParallel { SAM bam ->
                    
                    MeanCoverageEstimator meanEstimator = new MeanCoverageEstimator(bam, estRegions)
                    meanEstimator.sdThreshold = 0.5 // the default is slightly less accurate than wanted for this purpose
                    meanEstimator.minRegions = 30
                    double mean = meanEstimator.estimate()
                    
                    String sample = bam.samples[0]
                    log.info "Mean of $sample = $mean"
                    sampleMeans[sample] = mean
                }
                printer.setSampleMeans(sampleMeans)
            }
            else {
                log.info "Sample Means: " + printer.sampleMeans.collect { s, m -> "${s}=$m"}.join(", ")
            }
            log.info "Ordered Means: " + printer.orderedMeans.join(", ")
            
            bams.eachParallel { SAM bam ->
                combiner.processBAM(bam, this.scanRegions)
            } 
        }
            
        combiner << "stop"
        combiner.join()
            
        log.info "Stopping actor ..."
        printer << "stop"
        printer.join()
        log.info "Finished"
        
        if(phase1)
            return
        
        if(opts.stats) {
            log.info "Per sample coverage statistics are: " 
            Utils.table(out:System.err,
                ["Sample", "Mean (actual" + (opts.rel?",relative":"") + ")", "Mean (Est)", "StdDev"], 
                [ 
                    printer.samples, 
                    printer.samples.collect { printer.sampleStats[printer.samples.indexOf(it)].mean },
                    printer.samples.collect { printer.sampleMeans[it] },
                    printer.samples.collect { printer.sampleStats[printer.samples.indexOf(it)].standardDeviation }
                ].transpose()
            )
            if(opts.cv) {
                printCVInfo(printer.coeffvStats)
            }
        }
        
        if(opts.cvo) {
            new File(opts.cvo).withWriter { cvw ->
                writeCVOutput(cvw, printer.coeffvStats)
            }
        }
        
        if(opts.cvj) {
            new File(opts.cvj).withWriter { cvw ->
                writeCVOutputJS(cvw, printer.coeffvStats)
            }
        }
        
        if(opts.corr) {
            printCorrelationTable(printer)
        }
        
        if(opts.covo) {
            printCoverageJs(printer, opts.covo)
        }
        
        if(opts.gcprofile) {
            writeGCProfile(printer)
        }
    }
    
    private void printCoverageJs(CoveragePrinter printer, String fileName) {
        
        Map statsJson = [ 
            means: printer.samples.collectEntries { sample ->
                [ sample, printer.sampleStats[printer.samples.indexOf(sample)].mean] 
            },
                
            medians: printer.samples.collectEntries { sample ->
                [ sample, printer.rawCoverageStats[printer.samples.indexOf(sample)].median]
            },
        ]
        
        new File(fileName).withWriter { w ->
            w.println('covs = // NOJSON\n' + JsonOutput.prettyPrint(JsonOutput.toJson(statsJson)))
        }
        log.info "Wrote coverage stats to $fileName"
    }
    
    final List cvThresholds = (0..100).step(5)
    
    private void writeCVOutputJS(Writer w, IntegerStats coeffvStats) {
        w.println(
            'coeffvPercentiles = // NOJSON\n' + JsonOutput.toJson(
                [ 
                      cvThresholds.collect{ it/100},
                      cvThresholds.collect { thresh -> (1 - coeffvStats.fractionAbove(thresh)) }
                ].transpose()
            )
        )
    } 
    
    private void writeCVOutput(Writer w, IntegerStats coeffvStats) {
        w.println cvThresholds.collect{ it/100}*.toString().join('\t')
        w.println cvThresholds.collect { thresh -> Utils.perc(1 - coeffvStats.fractionAbove(thresh)) }.join('\t')
    }
    
    private void printCVInfo(IntegerStats coeffvStats) {
        log.info "Overall Coefficient of Variation: mean = ${coeffvStats.mean/100} median = ${coeffvStats.median/100}"
        Utils.table(out:System.err, topborder:true,
                    cvThresholds.collect{ it/100}*.toString(),
                    [cvThresholds.collect { thresh -> Utils.perc(1 - coeffvStats.fractionAbove(thresh)) }]
       )
    }

    private printCorrelationTable(CoveragePrinter printer) {
        
        final List<Stats> sampleStats = printer.sampleStats
        
        Matrix finalCorrelationMatrix = printer.correlationNumerators.transform { double value, int i, int j ->
            if(i==j)
                return 1.0
            else
                return value / (sampleStats[i].n * sampleStats[i].standardDeviation*sampleStats[j].standardDeviation)
        }
        
        finalCorrelationMatrix.@names = printer.samples
        
        println "Raw correlation matrix = \n"  + finalCorrelationMatrix
        
        finalCorrelationMatrix = new Matrix(finalCorrelationMatrix.getColumns(printer.correlationSamples))
        
        log.info "Computed correlations: "
       
        if(opts.cos) {
            opts.cos.each { String co ->
                log.info "Writing $co"
                new File(co).withWriter { w ->
                    if(co.endsWith('js')) 
                        writeCorrelationJs(w, printer, finalCorrelationMatrix)
                    else
                        writeCorrelationTSV(w, printer, finalCorrelationMatrix)
                }
            }
        }
        
        System.err.println "Final correlation: "
        finalCorrelationMatrix.@names = printer.samples
        finalCorrelationMatrix.sample = printer.correlationSamples
        System.err.println finalCorrelationMatrix.toString()
    }
    
    private writeGCProfile(CoveragePrinter printer) {
        
        log.info "Writing GC profile to $opts.gcprofile ..."
        
        List sampleGCProfiles = printer.gcBins.collect { sample, List<Stats> bins ->
            int binIndex = 0
            [
                sample: sample,

                gc: bins.collect {
                    Map info = [
                        bin: [binIndex*0.05, (binIndex+1)*0.05],

                        ratio: printer.gcBins[sample][binIndex].mean,

                        sd: printer.gcBins[sample][binIndex].standardDeviation
                    ]
                    ++binIndex
                    return info
                }
            ]
        }

        List interpolatedGC = sampleGCProfiles.collect { Map sampleGCInfo ->
            List<Map> gc = sampleGCInfo.gc
            return [ sample: sampleGCInfo.sample, gc: this.interpolateGCProfile(gc) ]
        }
        
        new File(opts.gcprofile).withWriter { w ->
            log.info "Computed GC profile: \n" + interpolatedGC
            w.println(JsonOutput.toJson(interpolatedGC))
        }            
        
    }
    
    List<Map> interpolateGCProfile(List<Map> gc) {
        
        log.info "Interpolating GC values: $gc"
        
        List<Map> usableBins = gc.grep {
            !it.ratio.isNaN()
        }
        
        log.info "Usable GC values are: $usableBins"
        
        // Treat the midpoint of each bin as the X coordinate for interpolation
        double [] x = [0d] + usableBins*.bin.collect { (it[0] + it[1]) / 2 } + [1.0d]
        double [] y = [0d] + usableBins*.ratio + [0d]
        
        LoessInterpolator interp = new LoessInterpolator()
        PolynomialSplineFunction gcFn = interp.interpolate(x, y)
        List<Map> result =  gc.collect { Map bin ->
            [
                bin: bin.bin,
                mean: Math.max(0,gcFn.value((bin.bin[0] + bin.bin[1])/2)),
                sd: bin.sd.isNaN() ? 0d : bin.sd
            ]
        }
        
        log.info "Interpolation result is $result"
        
        return result
    }

    private writeCorrelationTSV(BufferedWriter w, CoveragePrinter printer, Matrix finalCorrelationMatrix) {
        w.println printer.samples.join('\t')
        for(String sample in printer.correlationSamples) {
            w.println finalCorrelationMatrix[][printer.correlationSamples.indexOf(sample)].collect {fmt.format(it)}.join('\t')
        }
    }
    
    private writeCorrelationJs(BufferedWriter w, CoveragePrinter printer, Matrix finalCorrelationMatrix) {
        Map correlationBySample = printer.correlationSamples.collectEntries { corrSample ->
            [
                corrSample, 
                printer.samples.collectEntries { sample ->
                    [ 
                        sample, 
                        fmt.format(finalCorrelationMatrix[printer.samples.indexOf(sample)][printer.correlationSamples.indexOf(corrSample)]).toDouble()
                    ]
                }
            ]
        }
        w.println('corr=// NOJSON\n' + JsonOutput.prettyPrint(JsonOutput.toJson(correlationBySample)))
    }
     
}
