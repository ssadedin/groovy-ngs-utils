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
package gngs.coverage

import java.text.NumberFormat

import org.apache.commons.math3.stat.descriptive.SummaryStatistics

import gngs.*
import graxxia.IntegerStats
import graxxia.Matrix
import graxxia.Stats
import groovy.transform.CompileStatic
import groovy.util.logging.Log
import groovyx.gpars.GParsPool
import groovyx.gpars.actor.DefaultActor

@CompileStatic
class CoverageSummary {
    final PositionCounts countInfo
    final List<Double> values
    final Double coeffV
    final Double baseMean

    public CoverageSummary(final PositionCounts countInfo, final List<Double> values, final Double coeffV, final Double baseMean) {
        this.countInfo = countInfo;
        this.values = values;
        this.coeffV = coeffV;
        this.baseMean = baseMean
    }    
}

@Log
class CoverageSummarizer extends RegulatingActor<PositionCounts> {
    
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
    
    /**
     * Statistics for final computed coverage values - for example, these will be normalised
     * to the mean coverage etc. depending on the options supplied.
     */
    SummaryStatistics [] sampleStats = null
   
    List<String> sampleRegionStatsRegions = []
    
    Stats [] sampleRegionStats = null
    
    IntegerStats [] rawCoverageStats = null
    
    /**
     * Indexed by sample, a list of statistics for each target region
     */
    Map<String,List<Short>> sampleRegionMeans = null
    
    FASTA gcReference = null
    
    Region currentTarget = null
    
    double currentGc = -1.0d
    
    final NumberFormat numberFormat = NumberFormat.numberInstance
    
    final int numSamples
    
    boolean collectRegionStatistics = false
    
    boolean perBaseMean = false
    
    boolean outputBEDIds = false
    
    boolean outputTarget = false
    
    CoveragePrinter printer 
    
    CoverageSummarizer(Map options=[:], Writer w, List<String> samples) {
        super(50000,100000)
        this.progress = new ProgressCounter(withRate:true, log:log, withTime:true, extra: { "Computing stats on region: $currentTarget" })
        this.w = w
        this.samples = samples
        this.numSamples = samples.size()
        this.sampleStats = (1..samples.size()).collect { new SummaryStatistics() } 
        
        if(options.collectRegionStatistics) {
            this.collectRegionStatistics = true
            this.sampleRegionMeans = samples.collectEntries { [ it, new ArrayList(1000)] }
        }
        
        if(options.gcReference) {
            this.gcReference = options.gcReference
            log.info "Initializing GC profile bins"
            this.gcBins = this.samples.collectEntries {
                [ it, 
                  (1..20).collect { new Stats() }
                ]
            }
        }
        this.rawCoverageStats = (1..samples.size()).collect { new IntegerStats(1000) } 
        numberFormat.maximumFractionDigits=3
        numberFormat.minimumFractionDigits=0
        numberFormat.groupingUsed = false
        
        if(w != null) {
            this.printer = new CoveragePrinter(w)       
        
            if(options.outputTarget) {
                this.printer.outputTarget = true
            }
        }
//        this.printer.start()
        
//        this.batcher = new BatchedAcknowledgeableMessage(this.printer.pendingMessageCount, 20)
    }
    
    void setRelative(boolean value) {
        this.relative = value
        if(this.sampleMeans != null && !this.sampleMeans.isEmpty())
            this.setSampleMeans(this.sampleMeans)
    }
    
    void setSampleMeans(Map<String, Double> means) {
        
        this.sampleMeans = Collections.synchronizedMap(means)
        
        List missingSamples = this.samples.grep { !means.containsKey(it) }
        if(missingSamples) {
            throw new IllegalArgumentException("The following samples were missing from the provided set of means: " + missingSamples.join(',') +
                " (provided: ${means*.key.join(',')})")
        }
        
        // For efficiency, it is helpful to have the means pre-ordered in the same order 
        // that we want to output them in
        if(relative)
            orderedMeans = samples.collect { 1.0d } as double []
        else
            orderedMeans = samples.collect { sampleMeans[it] } as double []
    }
    
    Map<String, List<Stats>> gcBins 
    
    Map<String, Stats> currentGCBins
    
    int currentGCBin = -1
    
    @CompileStatic
    void process(PositionCounts countInfo) {
        final Region region = countInfo.region
        if(currentTarget != region) {
            newTargetRegion(region)
        }
        
        final List<Double> rawCounts = new ArrayList(samples.size())
        final Map<String,Double> counts = (Map<String,Double>)countInfo.counts
        for(int i=0; i<numSamples; ++i) {
            rawCounts[i] = counts[samples[i]]
        }
            
        List<Double> values 
        if(relative) {
            values = divideByMeans((Map<String,Double>)countInfo.counts)
            updateGCProfile(values)
        }
        else {
            values = rawCounts
        }
            
        updateRawCoverageStats(rawCounts)
        updateStats(values)
        
        if(std) {
            final Double valueMean = Stats.mean(values)
            values = values.collect { it /(0.01d +  valueMean) }
        }
        
        final Stats stats = (perBaseMean || coeffV) ? Stats.from(values) : null
        Double coeffVColumn = null
        if(coeffV) {
            final double coeffV = stats.standardDeviation / (1 + stats.mean)
            coeffVColumn = coeffV
            coeffvStats.addValue((int)(100*coeffV))
        }
        
        Double baseMean 
        if(perBaseMean) {
            baseMean = stats.mean
        }
        
//        writePosition(countInfo, values, coeffVColumn, baseMean)
        
        final String id = this.outputBEDIds ? this.currentTarget.extra : null
        
        if(!printer.is(null))
            printer.writePosition(countInfo, values, coeffVColumn, baseMean, id)
    }
    
    @CompileStatic
    List<Double> divideByMeans(final Map<String,Double> counts) {
       samples.collect{counts[it]/(1d + sampleMeans[it])} 
    }
    
    @CompileStatic
    void newTargetRegion(final Region region) {
        if(this.collectRegionStatistics) {
            updateRegionMeanCoverages()
            
            this.sampleRegionStats = new Stats[samples.size()].collect { new Stats() } as Stats[]
        }
        
        currentTarget = region
        if(gcReference) {
            currentGc = gcReference.gc(currentTarget)
            currentGCBin = (int)(currentGc * 100 / 5)
            // log.info "Calculated gc content $currentGc for region $currentTarget (bin $currentGCBin)"
        }
    }
    
    @CompileStatic
    void updateRegionMeanCoverages() {
        if(sampleRegionStats != null) {
            String targetValue = currentTarget?.toString()
            this.sampleRegionStatsRegions << targetValue
            [samples,sampleRegionStats].transpose().each { String sample, Stats stats -> 
                this.sampleRegionMeans[sample] << (Short)Math.round(stats.mean)
            }
        }
    }

    StringBuilder line = new StringBuilder(2048)
    
    BatchedAcknowledgeableMessage batcher 
    
    @CompileStatic
    void writePosition(final PositionCounts countInfo, final List<Double> values, final Double coeffV, final Double baseMean) {
//        CoverageSummary summary = new CoverageSummary(countInfo, values, coeffV)
//        this.printer.sendTo(summary)
//        batcher.batchTo(summary, printer)
        
//        printer.queue.offer(summary)
        
        if(w != null) {
            line.setLength(0)
            
            w.write(countInfo.chr)
            w.write('\t')
            w.write(countInfo.pos.toString())
            w.write('\t')
            if(coeffV != null) {
                w.write(numberFormat.format(coeffV))
                w.write('\t')
            }
            final int numValues = values.size()-1
            int i = 0;
            for(;i<numValues; ++i) {
                w.write(numberFormat.format(values[i]))
                w.write('\t')                                
            }
            w.write(numberFormat.format(values[numValues]))
            w.write('\n')
        }
    }
    
    @CompileStatic
    void updateRawCoverageStats(List<Double> values) {
        final int numValues = values.size()
        for(int i=0; i<numValues; ++i) {
            int intVal = (int)values[i]
            rawCoverageStats[i].addValue(intVal)
            if(collectRegionStatistics)
                sampleRegionStats[i].addValue(intVal)
        }
    } 
    
    @CompileStatic
    void updateGCProfile(List<Double> values) {
        if(currentGCBin>=0) {
            final int numSamples = values.size()
            for(int i=0; i<numSamples; ++i) {
                this.gcBins[samples[i]][currentGCBin].addValue(values[i])
            }
        }
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

    @Override
    public void onEnd() {
        log.info "Stopping summarizer at region $currentTarget"
        this.updateRegionMeanCoverages()
    }
    
    /**
     * Convenience method that computes coverage for the given alignments over the given regions,
     * and returns the resulting CoverageSummarizer.
     * <p>
     * Note: the base level coverage will be computed but ignored.
     * @param bams      
     * @param regions
     * @return
     */
    public static CoverageSummarizer summarize(List<SAM> bams, Regions regions) {
        List<String> samples = bams.collect { it.samples[0] }
        CoverageSummarizer summarizer = new CoverageSummarizer(null, samples, collectRegionStatistics: true)
        summarizer.start()
        CoverageCombinerActor combiner = new CoverageCombinerActor(summarizer, bams.size())
        combiner.start()
        
        GParsPool.withPool(bams.size()) { 
            bams.eachParallel { SAM bam ->
                CoverageCalculatorActor.processBAM(bam, regions, combiner, 1)
            }
        }
        log.info "Stopping ..."
        
        [combiner,summarizer].each { it.sendStop(); it.join() }
        
        return summarizer
    }
}

