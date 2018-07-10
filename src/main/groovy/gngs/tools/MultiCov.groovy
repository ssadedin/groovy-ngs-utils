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
import java.util.Iterator
import java.util.concurrent.atomic.AtomicInteger

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics

import gngs.*
import graxxia.IntegerStats
import graxxia.Matrix
import graxxia.Stats
import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import groovy.transform.ToString
import groovy.util.logging.Log
import groovyx.gpars.GParsPool
import groovyx.gpars.actor.Actor
import groovyx.gpars.actor.Actors
import groovyx.gpars.actor.DefaultActor
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SAMRecordIterator

@CompileStatic
@ToString
class SampleReadCount {
    String chr
    
    int pos 
    
    int reads
    
    String sample
}

@CompileStatic
class PositionReadCount {
    int position
    int reads
}

@Log
class CoverageReaderActor extends DefaultActor {
    
    AtomicInteger pending = new AtomicInteger(0)
    
    String currentChr
    
//    List<int[]> reads = new LinkedList()
    
    final OverlapTracker reads = new OverlapTracker()
    
    int pos
    
    Region currentRegion
    
    int currentReferenceIndex = -1
    
    Regions targetRegions
    
    Iterator<Region> regionIter
    
    Actor combiner
    
    String sample
    
    Set<Integer> bedReferenceIndexes
    
    List<String> bamContigs
    
    // for debug only
    SAMRecord currentRead
    
    ProgressCounter progress = new ProgressCounter(
        withRate:true, 
        extra:  {"region: $currentRegion, sample: $sample, readBuffer: ${reads.reads.size()}, pending: $pending"},
        log: log,
        timeInterval: 1000,
        lineInterval: 100
        )
    
    public CoverageReaderActor(SAM bam, Regions targetRegions, Actor combiner, String sample) {
        super();
        this.targetRegions = targetRegions
        this.regionIter = targetRegions.iterator();
        this.combiner = combiner;
        this.sample = sample;
        
        this.bamContigs = bam.getContigList()
        this.bedReferenceIndexes = targetRegions*.chr.unique().collect { String chr ->
            bamContigs.indexOf(chr)
        } as Set
        
        nextRegion()
    }

    @CompileStatic
    @Override
    void act() {
        loop { 
            react { msg ->
                if(msg == "stop") {
                    this.flushToEnd()
                    this.progress.end()
                    this.terminate()
                }
                else {
                    pending.decrementAndGet()
                    this.progress.count()
                    processRead((SAMRecord)msg)
                }
            }
        }
    }
    
    @CompileStatic
    void processRead(SAMRecord r) {
        
        
        if(r.referenceIndex < currentReferenceIndex)
            return
            
        if(currentRegion == null) // occurs if we ran out of regions => end
            return
        
        // Assumption: reads in BAM file 
        // are sorted the same as in the BED file
        this.currentRead = r
        assert r != null
        assert currentRegion != null
        
//        if(true)
//            return
//        
        while(r.referenceName != currentRegion.chr) {
            
            if(!bedReferenceIndexes.contains(r.referenceIndex)) {
                log.info "Ignoring contig $r.referenceName because not referenced in target regions"
                return
            }
            
            while(pos < currentRegion.to) {
                ++pos
                flushPosition()
            }
            if(regionIter.hasNext())
                if(!nextRegion())
                    return
        }
        
        flushToReadPosition(r)
        
        reads.add(r)
    }
    
    private void flushToEnd() {
        if(currentRegion == null)
            return
            
        while(currentRegion != null) {
            int regionEnd = currentRegion.to
            while(pos < regionEnd) {
                ++pos
                flushPosition()
            }
            
            nextRegion()
        }
    }
    
    @CompileStatic
    private flushToReadPosition(SAMRecord r) {
        final int alignmentStart = r.alignmentStart
        int regionEnd = currentRegion.to
        while(pos < alignmentStart) {
            ++pos
            if(pos > regionEnd) {
                nextRegion()
                
                // If switching to the next region caused us to hit a new 
                // chr then we have finished and can ignore this read
                // Note also that running out of regions causes a change in
                // currentReferenceIndex so that too will terminate here
                if(r.referenceIndex != currentReferenceIndex)
                    return
                    
                regionEnd = currentRegion.to
            }
            else
                flushPosition()
        }
    }
    
    @CompileStatic
    boolean nextRegion() {
       Region oldRegion = currentRegion
       
       if(!regionIter.hasNext()) {
           currentRegion = null
           currentReferenceIndex = -1
           return false
       }
           
       currentRegion = regionIter.next() 
       currentReferenceIndex = this.bamContigs.indexOf(currentRegion.chr)
       if(currentReferenceIndex < 0)
           throw new IllegalStateException("BED file sequence ${currentRegion.chr} is not found in BAM file sequences")
           
       pos = currentRegion.from
       if((oldRegion != null) && currentRegion.chr != oldRegion.chr) {
           log.info "Processing $currentRegion.chr"
           this.reads.reads.clear()
       }
       else {
           dropNonOverlapping()
       }
       return true
    }
    
    @CompileStatic
    void flushPosition() {
        dropNonOverlapping()
//        log.info "Flush $currentRegion.chr:$pos - ${reads.size()}"
        combiner << new SampleReadCount(chr: currentRegion.chr, pos: pos, reads: reads.reads.size(), sample: sample)
    }
    
    @CompileStatic
    void dropNonOverlapping() {
        this.reads.removeNonOverlaps(pos)
    }
    
    @CompileStatic
    void iteratorRemove(Iterator i) {
        i.remove()
    }
}

@Log
class CoverageCombinerActor extends DefaultActor {
    
    /**
     * How many samples to expect
     */
    int numSamples = 0
    
    String chr
    
    long pos = -1
    
    TreeMap<Long,Map<String,Integer>> counts = new TreeMap()
    
    AtomicInteger countSize = new AtomicInteger()
    
    Actor consumer
    
    ProgressCounter progress = new ProgressCounter(withRate: true, timeInterval: 1000, lineInterval: 200, log:log, 
        extra: { "Combining $chr:$pos (${XPos.parsePos(pos).startString()}) with ${counts[pos]}/$numSamples reported" })
    
    @Override
    void act() {
        loop {
            react { msg ->
                if(msg == "stop") {
                    this.progress.end()
                    log.info "Combiner terminating"
                    this.terminate()
                }
                else
                    processCount(msg)
            }
        }
    }
    
    @CompileStatic
    void processCount(SampleReadCount count) {
        
        Long xpos = XPos.computePos(count.chr, count.pos)
        if(pos == -1)
            pos = xpos
        
        Map<String,Integer> positionCounts = this.counts[xpos]
        if(positionCounts == null) {
            this.counts[xpos] = positionCounts = new HashMap()
        }
        
        positionCounts[count.sample] = count.reads
        progress.count()
        
        if(pos == xpos) {
            if(positionCounts.size() == numSamples) {
                        
                counts.pollFirstEntry()
                
                if(!counts.isEmpty())
                    pos = counts.firstKey()
                else
                    pos = -1
                consumer << [ chr: count.chr, pos: count.pos, counts: positionCounts]
            }
        }
    }
     
}

class CoveragePrinter extends DefaultActor {
    
    Writer w
    
    List<String> samples
    
    /**
     * Statistics for coefficient of variation - to avoid large memory use in tracking
     * individual values we (ab)use the integer stats class and bin coeffv as
     * integers from 0 - 100.
     */
    IntegerStats coeffvStats = new IntegerStats(100)
    
    /**
     * If set to true, the coverage value for each sample is divided by its mean
     */
    boolean relative
    
    /**
     * If true, the coverage at each base position will be divided by the mean of all
     * samples at the base position
     */
    boolean std
    
    /**
     * If set to true, the coefficient of variation of the coverage is printed as column 2
     */
    boolean coeffV
    
    /**
     * If set, correlation between this sample and other samples will be calculated and printed
     */
    List<String> correlationSamples = null
    
    /**
     * Cached indices of the samples to be used for correlation calculations
     */
    int [] correlationSampleIndices = null
    
    /**
     * Mean coverage for each sample
     * <p>
     * Only used if {@link #relative} is enabled
     */
    Map<String, Double> sampleMeans = Collections.synchronizedMap([:])
    
    Matrix correlationNumerators = null
    
    double [] orderedMeans = null
    
    Stats [] sampleStats = null
    
    final NumberFormat numberFormat = NumberFormat.numberInstance
    
    CoveragePrinter(Writer w, List<String> samples) {
        this.w = w
        this.samples = samples
        this.sampleStats = (1..samples.size()).collect { new Stats() } 
        numberFormat.maximumFractionDigits=3
        numberFormat.minimumFractionDigits=0
    }
    
    void act() {
        loop {
            react { msg ->
                if(msg == "stop")
                    terminate()
                else
                    writePosition(msg)
            }
        }
    }
    
    void setSampleMeans(Map<String, Double> means) {
        
        this.sampleMeans = means
        
        // For efficiency, it is helpful to have the means pre-ordered in the same order 
        // that we want to output them in
        if(relative)
            orderedMeans = samples.collect { 1.0d } as double []
        else
            orderedMeans = samples.collect { sampleMeans[it] } as double []
    }
    
    void writePosition(Map countInfo) {
        
        List<Double> values 
        if(relative) {
            values = samples.collect{countInfo.counts[it]/(1 + sampleMeans[it])}
        }
        else
            values = samples.collect{countInfo.counts[it]}
            
        updateStats(values)
        
        if(std) {
            Double valueMean = Stats.mean(values)
            values = values.collect { it /(0.01 +  valueMean) }
        }
        
        List coeffVColumn = []
        if(coeffV) {
            Stats stats = Stats.from(values)
            double coeffV = stats.standardDeviation / (1 + stats.mean)
            coeffVColumn = [numberFormat.format(coeffV)]
            coeffvStats.addValue((int)(100*coeffV))
        }
        
        w.println(([countInfo.chr, countInfo.pos] + coeffVColumn + values.collect{numberFormat.format(it)}).join('\t'))
    }
    
    @CompileStatic
    void updateStats(List<Double> values) {
        final int numValues = values.size()
        for(int i=0; i<numValues; ++i) {
            sampleStats[i].addValue(values[i])
        }
        computeCorrelation(values)
    }
    
    @CompileStatic
    void computeCorrelation(List<Double> values) {
        
        if(correlationSamples == null) 
            return 
            
        if(correlationSampleIndices == null) {
            correlationSampleIndices = correlationSamples.collect { samples.indexOf(it) } as int[]
            if(correlationSampleIndices.any { int i ->i <0 })
                throw new IllegalArgumentException("One or more samples $correlationSamples specified for correlation calculation is not found in samples for analysis: $samples")
        }
                
        if(this.correlationNumerators == null)
            this.correlationNumerators = new Matrix(new double[values.size()][values.size()])
            
        double [][] dataRef = this.correlationNumerators.matrix.dataRef
        
        for(int correlationSampleIndex in correlationSampleIndices) {
            double sampleDelta = values[correlationSampleIndex] - (relative ? 1.0d : sampleMeans[samples[correlationSampleIndex]])
            final int numValues = values.size()
            for(int i=0; i<numValues; ++i) {
                // TODO: we could only calculate one half of the matrix since it is symmetrical
                if(i != correlationSampleIndex) {
                    dataRef[i][correlationSampleIndex] += (values[i]  - orderedMeans[i]) * (sampleDelta)
                }
            }
        }
    }
}

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
            o 'Output file to write results to', longOpt: 'output', args:1
        }
    }
    
    @Override
    public void run() {
        
        fmt.maximumFractionDigits = 2
        
        resolveRegionsToScan()
        
        this.scanRegions = this.scanRegions.reduce().toSorted(new NumericRegionComparator())    
        
        this.bams = opts.arguments().collect { new SAM(it) }
        
        List<String> samples = bams*.samples*.getAt(0)
        log.info "Analysing coverage over ${Utils.humanBp(scanRegions.size())} for ${samples.join(',')}"
        
        def output = opts.o ? new File(opts.o) : System.out
        output.withWriter { w ->
            
            CoveragePrinter printer = new CoveragePrinter(w, samples)
            if(opts.cv)
                printer.coeffV = true
            if(opts.rel)
                printer.relative = true
            if(opts.std)
                printer.std = true
                
            if(opts.corr)
                printer.correlationSamples = opts.corr == '.' ? samples : opts.corr.tokenize(',')
                
            run(w, printer)
        }
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
    
    void run(Writer w, CoveragePrinter printer) {
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
            
            bams.eachParallel { SAM bam ->
                processBAM(bam, combiner)
            } 
        }
            
        combiner << "stop"
        combiner.join()
            
        log.info "Stopping actor ..."
        printer << "stop"
        printer.join()
        log.info "Finished"
        
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
                        fmt.format(finalCorrelationMatrix[printer.samples.indexOf(sample)][printer.correlationSamples.indexOf(corrSample)])
                    ]
                }
            ]
        }
        w.println('corr=// NOJSON' + JsonOutput.prettyPrint(JsonOutput.toJson(correlationBySample)))
    }
     
    /**
     * Last time a warning was issued about throttled reading - only want to log this
     * once per minute or so
     */
    long throttleWarningMs = 0
    
    @CompileStatic
    void processBAM(SAM bam, CoverageCombinerActor combiner) {
        
        String sample = bam.samples[0]
        
        CoverageReaderActor cra = new CoverageReaderActor(bam, scanRegions, combiner, sample) 
        cra.start()
        
        List<String> chrs = this.scanRegions*.chr.unique()
        for(String chr in chrs) {
            int start = Math.max(0, this.scanRegions.index[chr].ranges.firstKey() - 1000)
            int end = this.scanRegions.index[chr].ranges.lastKey() + 1000
            log.info "Scan $chr from $start to $end"
            bam.withIterator(new Region(chr, start, end))  { SAMRecordIterator iter ->
                while(iter.hasNext()) {
                    cra << iter.next()
                    
                    int pending = cra.pending.incrementAndGet()
                    if(pending > 50000) {
                        long nowMs = System.currentTimeMillis()
                        if(nowMs - throttleWarningMs > 60000) {
                            log.info "Throttling reading due to downstream congestion"
                            throttleWarningMs = nowMs
                        }
                        Thread.sleep(50)
                    }
                }
            }
        }
        log.info "Sending stop message to CRA ${bam.samples[0]} ..."
        cra << "stop"
        cra.join()
    }
}
