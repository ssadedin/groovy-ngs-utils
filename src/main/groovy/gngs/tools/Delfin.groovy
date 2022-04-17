package gngs.tools

import java.util.concurrent.ConcurrentHashMap
import java.util.logging.Level

import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction

import com.twosigma.beakerx.chart.histogram.Histogram

import gngs.*
import graxxia.*
import groovy.transform.CompileStatic
import groovy.util.logging.Log
import groovyx.gpars.GParsPool
import groovyx.gpars.GParsPoolUtil
import groovyx.gpars.actor.Actor
import groovyx.gpars.actor.Actors

/**
 * A simple exon-by-exon deletion finder
 * 
 * @author Simon Sadedin
 */
@Log
class Delfin extends ToolBase {
    
    final static int MINIMUM_TARGET_REGIONS = 20
    
    /**
     * Likelihood ratio at which to output deletions
     */
    double deletionCallThreshold = 3.5

    /**
     * Likelihood ratio at which to output duplications
     */
    double duplicationCallThreshold = 3.5
    
    /**
     * The maximum number of PCA components that will be removed
     */
    int maxPCAComponents = 20
    
    /**
     * An indicator of the number deletion calls expected per target region
     * 
     * (or, approximately, the probability any given target is deleted state, though
     * that interpretation does not account for multi-target calls)
     */
    double deletionRate = 30 / 7088
    
    /**
     * The expected number of deletion calls is calculated from the per-target region rate
     * but this can lead to an unreasonably low bar for small chromosomes
     */
    int minExpectedDeletionCallsPerChromosome = 5
    
    Regions targets
    
    Matrix covs
    
    Map<String,Regions> results = new ConcurrentHashMap()
    
    Writer optimisationWriter
    
    /**
     * row indices of the samples being tested for CNV calls in the matrix 
     */
    private List<Integer> testSampleIndices = null

    /**
     * row indices of the samples being used only as controls for CNV calls in the matrix 
     */
    private List<Integer> controlSampleIndices = null
    
    private String debugChr = null
    
    private int debugRegionIndex = -1
    
    private FASTA reference
    
    /**
     * Matrix containing a row for each sample and a column for each GC bin,
     * with the relative inflation or deflation of coverage for that region
     * based on GC content, for the given sample.
     */
    private Matrix gcNormalisationFactors = null

    static void main(args) {
        cli('-t <target regions bed> -cov <cov1> [-cov <cov2> ...] -o <deletions.tsv>', args) {
            t 'Target region bed file', args:1, required: true, type: File, longOpt: 'targets'
            o 'Output TSV file containing deletions (stdout)', args:1, required: false, type: File
            prior 'Prior probability that target region is deleted (0.004)', args: 1, required: false, type: Double
            lrt 'Likelihood ratio threshold to emit CNV calls at (3.5)', args:1, required: false, type: Double
            del_lrt 'Likelihood ratio threshold to emit deletion calls at (3.5)', args:1, required: false, type: Double
            dup_lrt 'Likelihood ratio threshold to emit duplication calls at (3.5)', args:1, required: false, type: Double
            maxpc 'Set maximum principle components to remove from data (20)', args:1, required: false, type: Integer
            lr 'Output file to save likelihood ratio estimates in', args:1, required: false, type: File
            ref 'Reference fasta; if provided, GC normalisation will be used', args:1, required: false, type: String
            s 'Samples to test for CNVs (all)', args:'*', required: false, type:String
            chr 'Process the given chromosome only (false, process all)', args: 1, required: false
            cov 'Coverage depth file, in GATK format, tab separated with a column of read count for each target region', args: '*', type: File, required: true
            optout 'File to save optimisation metrics to (disabled)', args:1, required: false, type: File
            dr 'Debug region: output verbose statistics for given target region', args:1, required: false
        }
    }

    @Override
    public void run() {
        
        log.info "Loading $opts.t"
        targets = new BED(opts.t).load(withExtra:true)
        
        log.info "Loaded ${targets.numberOfRanges} regions (${Utils.humanBp(targets.size())} from $opts.t"
        
        if(opts.optout)
            this.optimisationWriter = opts.optout.newWriter()

        initialiseDebugRegion()
       
        GParsPool.withPool(6) {
            List<Matrix> covList = opts.covs.collectParallel { loadCoverage(it) }
            covs = Matrix.concat(covList)
        }
        
        log.info "Loaded ${covs.size()} coverage data sets"
        
        if(opts.ref) {
            log.info "Enabling GC normalisation using reference $opts.ref"
            this.reference = new FASTA(opts.ref)
            initGCNormalisation()
        }
        
        List<String> chrs = opts.chr ? [opts.chr] : targets.allRanges*.key
        
        log.info "Processing ${chrs.size()} chromosomes separately in parallel"
        
        if(opts.lrt != false) {
            this.deletionCallThreshold = opts.lrt
            this.duplicationCallThreshold = opts.lrt
        }
        
        if(opts.del_lrt != false) {
            this.deletionCallThreshold = opts.del_lrt
        }
        
        if(opts.dup_lrt != false) {
            this.duplicationCallThreshold = opts.dup_lrt
        }
        
        if(opts.maxpc) {
            this.maxPCAComponents = opts.maxpc
        }
        
        if(opts.prior != false)
            this.deletionRate = opts.prior
        
        this.testSampleIndices = covs.sample.findIndexValues { it in opts.ss }

        this.controlSampleIndices = covs.sample.findIndexValues { !(it in opts.ss) }
        
        if(!controlSampleIndices) {
            log.info "There are no control samples so test samples are being used as controls"
            this.controlSampleIndices = this.testSampleIndices
        }
        
        List<Exception> errors = []

        List actors = chrs.collect { chr ->
            Actors.actor {  
                
                delegate.metaClass.onException = { e ->
                    log.log Level.SEVERE, "Error occurred processing $chr: ", e
                    errors << e
                }
                
                Regions chrResults = analyseChromosome(chr)
                results[chr] = chrResults
            } 
        }
        
        log.info "Started actors to process each chromosome"
        
        actors*.join()
        
        if(errors) {
            log.severe "One or more errors occurred during analysis - exiting with error status"
            System.exit(1)
        }
        
        saveResults()
        
        if(this.optimisationWriter)
            this.optimisationWriter.close()
        
        log.info "Finished"
    }

    Region debugRegion 

    private initialiseDebugRegion() {
        if(opts.dr) {
            this.debugRegion = new Region(opts.dr)
            debugRegionIndex = (targets.findAll { it.chr == debugRegion.chr }).findIndexOf { it.overlaps(debugRegion) }
            if(debugRegionIndex<0)
                throw new IllegalArgumentException("No overlap of debug region $debugRegion could be found with target regions in $opts.t")
            debugChr = debugRegion.chr
            log.info "Debug region $debugRegion identified at index $debugRegion in $debugChr"
        }
    }

    private void saveResults() {
        log.info "Saving results to $opts.o"
        Utils.writer(opts.o).withWriter { w ->
            w.write(['chr','start','end','sample','lr', 'type'].join('\t'))
            w.write('\n')
            results.each { chr, cnvs ->
                cnvs.each { Region cnv ->
                    w.write([chr,cnv.from,cnv.to,cnv.sample,cnv.lr, cnv.type].join('\t'))
                    w.write('\n')
                }
            }
        }
    }
    
    final Regions analyseChromosome(final String chr) {
        
        List<Integer> colIndexes = targets.findIndexValues { it.chr == chr }

        log.info "Analysing $chr (${colIndexes.size()} target regions)"
        
        assert colIndexes.size() >= MINIMUM_TARGET_REGIONS && colIndexes.size() > covs.rowDimension : 
            "Chromosome $chr has only ${colIndexes.size()} target regions which is too few to analyse correctly for deletions"
            
        Matrix subset = covs[][colIndexes] 
        subset.@displayColumns = 10
        subset.@displayRows = 8
        subset.@displayPrecision = 2
        
        assert subset != null
        
        log.info "Coverages for $chr are:\n$subset"
        
        Matrix std = subset.standardise()
        
        log.info "Std for $chr are:\n$std"
        
        log.info "Calculating PCA basis vectors on $chr ..."
        
        if(maxPCAComponents>subset.rowDimension) {
            log.info "Too few samples to apply max pca componetns $maxPCAComponents: will reduce to ${subset.rowDimension-1}"
            maxPCAComponents = subset.rowDimension-1
        }
        
        Matrix reduced = std.reduce(maxPCAComponents)

        // Note: important that samples listed in order they are found in cov matrix
        List<String> samples = opts.ss ? covs.sample.findAll { it in opts.ss }: covs.sample 
        
        if(!samples)
            throw new IllegalStateException("None of the specified samples in $opts.ss where found in the provided coverage depth tables: \n" + covs.sample.join('\n') )
        
        List<List<List<Double>>> delAndDupCallsBySample = samples.collect { 
            callAndOptimise(chr, it, subset, std, reduced) 
        }

        Matrix deletionLRs = delAndDupCallsBySample.collect { 
            it[0]
        } as Matrix
        deletionLRs.@names = subset.@names
        deletionLRs.sample = samples
        
        Matrix duplicationLRs = delAndDupCallsBySample.collect { 
            it[1]
        } as Matrix
        duplicationLRs.@names = subset.@names
        duplicationLRs.sample = samples
        
        
        def allDels = deletionLRs*.findIndexValues { it > deletionCallThreshold }
        def allDups = duplicationLRs*.findIndexValues { it > duplicationCallThreshold }
        
        Regions chrRegions = targets[colIndexes] as Regions

        Regions mergedDels = segmentCNVs(samples, allDels, chr, deletionLRs, chrRegions)
        mergedDels.each { it.type = 'DEL' }

        Regions mergedDups = segmentCNVs(samples, allDups, chr, duplicationLRs, chrRegions)
        mergedDups.each { it.type = 'DUP' }

        return (mergedDels + mergedDups) as Regions
    }

    /**
     * Process the given pair of lists containing (samples, corresponding indexes of target regions called) and
     * (a) convert the target region indexes to Region objects, and (b) identify 
     * regions where multiple calls are adjacent and merge them to a single Region.
     */
    private Regions segmentCNVs(List samples, List<List> allCalls, String chr, Matrix allLRs, Regions chrRegions) {
        Regions mergedCNVs =
                [samples, allCalls].transpose().collectMany { String sampleId, List<Integer> regionIndices ->
                    log.info "Found ${regionIndices.size()} for $sampleId in $chr"
                    int sampleIndex = allLRs.sample.indexOf(sampleId)

                    Regions sampleUnmergedCNVs = chrRegions[regionIndices].eachWithIndex { cnv, index ->
                        cnv.sample = sampleId
                        cnv.lr = allLRs[sampleIndex][(int)regionIndices[index]]
                        cnv.index = regionIndices[index]
                    }

                    // Merge consecutive ranges
                    return mergeConsecutiveRegions(sampleUnmergedCNVs)
                }
        return mergedCNVs
    }
    
    /**
     * Tries calling deletions while removing increasing numbers of prinicipal components. When
     * an acceptably small number is returned, stops and returns both duplications and deletions.
     * 
     * @return List containing two elements: [indexes of called deletions, indexes of called duplications]
     */
    List<List<Double>> callAndOptimise(String chr, String sample, Matrix scaled, Matrix std,  Matrix reduced) {
        
        log.info "Calculating likelihoods for $sample on $chr"

        // The index of the sample to call CNVS for in the main matrix
        int sampleIndex = covs.sample.findIndexOf { it == sample }
        
        // Deletion rate is expressed as prior probability that any specific region
        // is deleted
        int maxDeletionCalls = std.columnDimension * deletionRate

        log.info "Expect $maxDeletionCalls calls based on deletionRate=$deletionRate and ${std.columnDimension} target regions in $chr"
        
        if(maxDeletionCalls < minExpectedDeletionCallsPerChromosome) {
            log.info "Expected number of deletion calls adjusted up from $maxDeletionCalls to $minExpectedDeletionCallsPerChromosome to meet minimum threshold"
            maxDeletionCalls = minExpectedDeletionCallsPerChromosome
        }

        List deletionLRs 
        List dupLRs 
        int nComponentsToRemove = 0
        int cnvCount = 0
        int minCNVCount = 0
        List<Matrix> cnvLRs
        List<Matrix> bestCnvLRs
        int bestCount = Integer.MAX_VALUE
        while(!deletionLRs || cnvCount > maxDeletionCalls) {

            cnvLRs = normaliseAndComputeLikelihoods(sample, chr, scaled, std, reduced, nComponentsToRemove)
            
            Matrix deletionMatrix = cnvLRs[0]
            
            deletionLRs = deletionMatrix[sampleIndex] as List
            dupLRs = cnvLRs[1][sampleIndex] as List

            nComponentsToRemove++
            
            int delCount = deletionLRs.count { it > deletionCallThreshold}
            int dupCount = dupLRs.count { it > duplicationCallThreshold }
            
            cnvCount = delCount  + dupCount
            
            if(Math.abs(cnvCount - maxDeletionCalls) < (Math.abs(bestCount - maxDeletionCalls))) {
                bestCount = cnvCount
                bestCnvLRs = cnvLRs
            }
           
            log.info "Analysis for $sample/$chr with $nComponentsToRemove components removed produced $cnvCount calls ($delCount dels, $dupCount dups)"
            
            if(this.optimisationWriter)
                this.optimisationWriter.write("$cnvCount\t$maxDeletionCalls\t$nComponentsToRemove\t$deletionCallThreshold")

            if(nComponentsToRemove>=maxPCAComponents) {
                log.info "Exceeded maximum $maxPCAComponents PCA components removed - using optimum result producing $bestCount calls (targeting $maxDeletionCalls calls)"
                if(bestCnvLRs != null) {
                    cnvLRs = bestCnvLRs
                }
                deletionLRs = cnvLRs[0][sampleIndex] as List
                dupLRs = cnvLRs[1][sampleIndex] as List
                break
            }
        }

//        if(opts.lr) {
//            Utils.writer("${opts.lr}.${sample}.${chr}.tsv").withWriter { w -> 
//                sampleLRs.save(w, r:true) 
//            }
//        }

        // Identify candidate CNV regions
        def deletionCalls = deletionLRs.findIndexValues { it > deletionCallThreshold }
        
        log.info "Deletion lrs are: \n$deletionCalls"
       

        def dupCalls = dupLRs.findIndexValues { it > duplicationCallThreshold }

        return [deletionLRs,dupLRs]
    }
 
    /**
     * Transform the matrix of row-normalised coverage values into a matrix 
     * of log-scaled likelihood ratios comparing the z-score of a diploid
     * model to a single-ploidy model.
     * 
     * @param row-normalised coverage values, row per sample, column per target region
     * @return  matrix of likelihoods
     */
    @CompileStatic
    List<Matrix> normaliseAndComputeLikelihoods(final String sample, final String chr, final Matrix scaled, final Matrix std, final Matrix reduced, final int nComponentsToRemove) {

        assert std != null

       
        Matrix basis = (Matrix)reduced.metadata.basis
        
        log.info "Calculating likelihoods with $nComponentsToRemove components removed for $sample $chr"
        Matrix approxStd
        Matrix approx
        if(nComponentsToRemove > 0) {
            approx = std - (Matrix)((reduced.columns as List)[0..<nComponentsToRemove] as Matrix) * (Matrix)(basis[0..<nComponentsToRemove])
            approxStd = approx.standardise()
        }
        else {
            approx = std
            approxStd = std
        }
        
        if(debugChr == chr) {
            dumpDebug("Standardised values", approxStd)
        }
        
        log.info "Calculate deletion likelihoods for $sample (pcomp=$nComponentsToRemove)..."
        Matrix deletionLRs = computeLikelihoods(scaled, approxStd, std, sample, chr, 0.5d)

        log.info "Calculate duplication likelihoods for $sample (pcomp=$nComponentsToRemove)..."
        Matrix duplicationLRs = computeLikelihoods(scaled, approxStd, std, sample, chr, 1.5d)

        return [deletionLRs,duplicationLRs]
    }

    @CompileStatic
    private Matrix computeLikelihoods(Matrix scaled, Matrix approxStd, Matrix std, String sample, String chr, final double cnExpectedValue) {

        // Diploid regions are by definition std dev 1 because the matrix was standardised
        FastNormal diploid = new FastNormal(0, 1)

        // For deletions we must calculate a separate estimated mean from each column because the
        // standard deviation of each column is different
        List<FastNormal> cnvDistributions = scaled.columns.collect {  MatrixColumn c ->
            Stats colStats = Stats.from((List<Double>)c[controlSampleIndices])
            if(colStats.standardDeviation>0) {
                new FastNormal((cnExpectedValue * colStats.mean - colStats.mean) / colStats.standardDeviation, 1)
            }
            else {
                // Most likely all zeros in the control samples
                null
            }
        }

        Matrix lrs = approxStd.transform { double x, int i, int j -> if(cnvDistributions[j]!=null) cnvDistributions[j].logDensity(x) else Double.MIN_VALUE } -
        approxStd.transform { double x -> diploid.logDensity(x) }

        lrs.@names = std.@names

        if(debugChr == chr) {
            logStats(sample, chr, scaled, approxStd, cnvDistributions, diploid, cnExpectedValue)
        }
        return lrs
    }
    

    private logStats(String sample, String chr, Matrix approx, Matrix approxStd, List<FastNormal> cnvPredictedDistributions, FastNormal diploid, double cnExpectedValue) {
        
        FastNormal cnvDist = cnvPredictedDistributions[debugRegionIndex]
        log.info "CNV (mean, stddev) for $sample in $chr at $debugRegionIndex (predicted for CNV: $cnExpectedValue): " + [
            cnvDist.mean,
            cnvDist.standardDeviation
        ].join(', ')
        
        double stdValue = approxStd[((List<String>)approxStd.getProperty('sample')).indexOf(sample)][debugRegionIndex]
        log.info "Standardised value of $sample at $debugRegion is " + stdValue
            
        double diploidDensity = diploid.logDensity(stdValue)
        log.info "Diploid density = " + diploidDensity
        
        double cnvDensity = cnvDist.logDensity(stdValue)

        log.info "Sample density (likelihood) given CNV = " + cnvDensity
        log.info "Likelihood ratio = " + (cnvDensity - diploidDensity)
        if((cnvDensity - diploidDensity) > this.deletionCallThreshold) {
            log.info "CNV call WILL be emitted for $opts.dr"
        }
        else {
            log.info "CNV call NOT emitted for $opts.dr"
        }
    }
    
    private dumpDebug(msg, Matrix approxStd) {
        IntRange debugRange = (debugRegionIndex-5)..(debugRegionIndex+5)
        Matrix subset = approxStd[][debugRange]
        subset.@displayRows = 100
        subset.sample = approxStd.sample
        log.info("$msg at debug region $debugRegionIndex: \n" + subset)

    }

    private List<Region> mergeConsecutiveRegions(final Regions unmergedCNVs) {
        return unmergedCNVs.inject([]) { List<Region> result, Region cnv  ->
            if(result && cnv.index == (result[-1].index+1)) {
                // Merge this CNV with the adjacent call
                Region adjacent = result[-1]
                
                log.info "Merge $adjacent ($adjacent.index) with $cnv ($cnv.index)"
                Region merged = new Region(cnv.chr, adjacent.from, cnv.to)
                merged.sample = cnv.sample
                merged.lr = Math.max(adjacent.lr, cnv.lr )
                merged.index = cnv.index
                result[-1] = merged
            }
            else
                result << cnv
            result
        }        
    }
    
    /**
     * Loads per-target-region coverage depth from the given provided files, which 
     * are in the form of an initial column containing the sample and then 
     * a column for each target region containing the read depth over that exon.
     * 
     * @param path  path to the bed file
     * @return  a Matrix with a row per sample, normalised so that the mean of each row is 1
     */
    Matrix loadCoverage(File path) {
        
        log.info "Loading coverage data from $path"
        
        Utils.reader(path) { r ->
        
            def controls = Matrix.load(r)
                             .normaliseRows()
                             .fillna(0)
                             
            assert controls.columnDimension == this.targets.numberOfRanges : 
                "Provided coverage file $path has different number of target regions ($controls.columnDimension) to the provided bed file $opts.t ($targets.numberOfRanges)"
            
            if(controls.properties.containsKey('Mean')) { // gngs bug
                controls.sample = controls.Mean
                controls.properties.remove('Mean')
            }
            controls.names = (0..<controls.columnDimension).collect { "T" + it}
            controls.@displayColumns = 12
            return controls
        }
    }
    
    
    double [] targetRegionGCRfractions
    
    @CompileStatic
    void initGCNormalisation() {
        log.info "Calculating GC content of target regions ..."
        this.targetRegionGCRfractions = new double[this.targets.numberOfRanges]

        this.calculateIntervalGCs()
        
        log.info "Finished calculating ${this.targets.numberOfRanges} GC fractions"

        final Matrix normCovs = this.covs.normaliseRows()
        
        log.info "Calculating GC factors for each target region ..."
        Matrix allRelGCs = normCovs.transformRows { double [] values, int rowIndex ->
            Binner binner = new Binner(20, 0d, 1d)
            
            log.info "Calculating GC curve for ${values.size()} regions for sample $rowIndex (on ${targetRegionGCRfractions.size()} target regions)"

            List<Stats> stats = binner.stats(this.targetRegionGCRfractions, values)

            PolynomialSplineFunction fn = binner.splineFromStats([startPoint:[0d,0d], endPoint:[1.0d,0d]], stats)
            
            final double [] relGCs = binner.midPoints.collect { fn.value(it) } as double[]
            return relGCs
        }
        
       
        log.info "Finished calculating GC normalisation factors: \n" + allRelGCs.toString()

        this.gcNormalisationFactors= allRelGCs.normaliseColumns().transform { double x -> Double.isNaN(x) ? 0d : x }
        
        this.gcNormalisationFactors.setProperty('sample',this.covs.getProperty('sample'))

        log.info "Relative adjustment factors for GC: \n" + gcNormalisationFactors.toString()
     }
    
    void calculateIntervalGCs() {
        def targetsFile = opts.t
        def gcCacheFile = new File("${targetsFile.name}.gc.tsv") 
        if(gcCacheFile.exists()) {
            log.info "Using cached GC regions for target file $targetsFile"
            this.targetRegionGCRfractions = gcCacheFile.readLines()*.toDouble() as double[]
            return
        }

        GParsPool.withPool(2) {
            use(GParsPoolUtil) {
                this.targets.eachWithIndexParallel { target, i ->  
                    this.targetRegionGCRfractions[i] = this.reference.gc(target)
                } 
            }
        }
        
        gcCacheFile.withWriter { w ->
            for(gc in this.targetRegionGCRfractions) {
                w.write(String.valueOf(gc))
                w.write('\n')
            }
        }
    }
}
