package gngs.tools

import java.util.concurrent.ConcurrentHashMap

import gngs.*
import graxxia.*
import groovy.transform.CompileStatic
import groovy.util.logging.Log
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
     * The maximum number of PCA components that will be removed
     */
    int maxPCAComponents = 20
    
    /**
     * Treat more than this number of deletion calls within a chromosome as a sign that
     * batch effects are creating artefactual CNVs
     */
    int maxDeletionCallsPerChromosome = 30
    
    Regions targets
    
    Matrix covs
    
    Map<String,Regions> results = new ConcurrentHashMap()
    
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

    static void main(args) {
        cli('-t <target regions bed> -cov <cov1> [-cov <cov2> ...] -o <deletions.tsv>', args) {
            t 'Target region bed file', args:1, required: true, type: File, longOpt: 'targets'
            o 'Output TSV file containing deletions (stdout)', args:1, required: false, type: File
            lrt 'Likelihood ratio threshold to emit deletion calls at (3.5)', args:1, required: false, type: Double
            maxpc 'Set maximum principle components to remove from data (20)', args:1, required: false, type: Integer
            lr 'Output file to save likelihood ratio estimates in', args:1, required: false, type: File
            s 'Samples to test for CNVs (all)', args:'*', required: false, type:String
            chr 'Process the given chromosome only (false, process all)', args: 1, required: false
            cov 'Coverage depth file, in GATK format, tab separated with a column of read count for each target region', args: '*', type: File, required: true
            dr 'Debug region: output verbose statistics for given target region', args:1, required: false
        }
    }

    @Override
    public void run() {
        
        log.info "Loading $opts.t"
        targets = new BED(opts.t).load(withExtra:true)
        
        log.info "Loaded ${targets.numberOfRanges} regions (${Utils.humanBp(targets.size())} from $opts.t"
        
        initialiseDebugRegion()
        
        covs = Matrix.concat(*opts.covs.collect { loadCoverage(it) })
        
        log.info "Loaded ${covs.size()} coverage data sets"
        
        List<String> chrs = opts.chr ? [opts.chr] : targets.allRanges*.key
        
        log.info "Processing ${chrs.size()} chromosomes separately in parallel"
        
        if(opts.lrt != false) {
            this.deletionCallThreshold = opts.lrt
        }
        
        if(opts.maxpc) {
            this.maxPCAComponents = opts.maxpc
        }
        
        this.testSampleIndices = covs.sample.findIndexValues { it in opts.ss }

        this.controlSampleIndices = covs.sample.findIndexValues { !(it in opts.ss) }
        
        if(!controlSampleIndices) {
            log.info "There are no control samples so test samples are being used as controls"
            this.controlSampleIndices = this.testSampleIndices
        }

        List actors = chrs.collect { chr ->
            Actors.actor {  
                Regions chrResults = analyseChromosome(chr)
                results[chr] = chrResults
            } 
        }
        
        log.info "Started actors to process each chromosome"
        
        actors*.join()
        
        saveResults()
        
        log.info "Finished"
    }

    private initialiseDebugRegion() {
        if(opts.dr) {
            Region debugRegion = new Region(opts.dr)
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
            w.write(['chr','start','end','sample','lr'].join('\t'))
            w.write('\n')
            results.each { chr, cnvs ->
                cnvs.each { Region cnv ->
                    w.write([chr,cnv.from,cnv.to,cnv.sample,cnv.lr].join('\t'))
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
        Matrix reduced = std.reduce(maxPCAComponents)
       
        // Note: important that samples listed in order they are found in cov matrix
        List<String> samples = opts.ss ? covs.sample.findAll { it in opts.ss }: covs.sample 
        
        Matrix lrs = samples.collect { callAndOptimise(chr, it, subset, std, reduced) } as Matrix
        lrs.@names = subset.@names
        lrs.sample = samples
        
        def allCNVs = lrs*.findIndexValues { it > deletionCallThreshold }
        
        Regions chrRegions = targets[colIndexes] as Regions
        Regions mergedCNVs =
            [samples, allCNVs].transpose().collectMany { String sampleId, List<Integer> regionIndices ->
                log.info "Found ${regionIndices.size()} for $sampleId in $chr"
                int sampleIndex = lrs.sample.indexOf(sampleId)
                
                Regions sampleUnmergedCNVs = chrRegions[regionIndices].eachWithIndex { cnv, index ->
                    cnv.sample = sampleId 
                    cnv.lr = lrs[sampleIndex][(int)regionIndices[index]]
                    cnv.index = regionIndices[index]
                }
                
                // Merge consecutive ranges
                return mergeConsecutiveRegions(sampleUnmergedCNVs)
            }
         

        return mergedCNVs
    }
    
    List<Double> callAndOptimise(String chr, String sample, Matrix scaled, Matrix std,  Matrix reduced) {
        
        log.info "Calculating likelihoods for $sample on $chr"

        // The index of the sample to call CNVS for in the main matrix
        int sampleIndex = covs.sample.findIndexOf { it == sample }

        List lrs 
        int nComponentsToRemove = 0
        int cnvCount = 0
        while(!lrs || cnvCount > maxDeletionCallsPerChromosome) {
            lrs = computeLikelihoods(sample, chr, scaled, std, reduced, nComponentsToRemove)[sampleIndex] as List
            nComponentsToRemove++
            
            cnvCount = lrs.count { it > deletionCallThreshold} 
            
            log.info "Analysis for $sample/$chr with $nComponentsToRemove components removed produced $cnvCount calls"

            if(nComponentsToRemove>=maxPCAComponents)
                break
        }

//        if(opts.lr) {
//            Utils.writer("${opts.lr}.${sample}.${chr}.tsv").withWriter { w -> 
//                sampleLRs.save(w, r:true) 
//            }
//        }

        // Identify candidate CNV regions
        def cnvs = lrs.findIndexValues { it > deletionCallThreshold }
        
        log.info "CNV lrs are: \n$cnvs"
       
        return lrs
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
    Matrix computeLikelihoods(final String sample, final String chr, final Matrix scaled, final Matrix std, final Matrix reduced, final int nComponentsToRemove) {

        assert std != null

       
        Matrix basis = (Matrix)reduced.metadata.basis
        
        log.info "Calculating likelihoods with $nComponentsToRemove components removed for $sample $chr"
        Matrix approxStd
        if(nComponentsToRemove > 0) {
            Matrix approx = std - (Matrix)((reduced.columns as List)[0..<nComponentsToRemove] as Matrix) * (Matrix)(basis[0..<nComponentsToRemove])
            approxStd = approx.standardise()
        }
        else {
            approxStd = std
        }
        
        if(debugChr == chr) {
            dumpDebug("Standardised values", approxStd)
        }
        
        // Diploid regions are by definition std dev 1 because the matrix was standardised
        FastNormal diploid = new FastNormal(0, 1)

         // For deletions we must calculate a separate estimated mean from each column because the 
        // standard deviation of each column is different
        List<FastNormal> deletions = scaled.columns.collect {  MatrixColumn c ->
            Stats colStats = Stats.from((List<Double>)c[controlSampleIndices])
            if(colStats.standardDeviation<0) {
                log.info "Wot? $colStats.standardDeviation"
            }

            new FastNormal(-0.5 * colStats.mean / colStats.standardDeviation, colStats.standardDeviation) 
        }
        
        if(debugChr == chr) {
            logStats(sample, chr, approxStd, deletions, diploid)
        }
        
        Matrix lrs = approxStd.transform { double x, int i, int j -> deletions[j].logDensity(x) } - 
                     approxStd.transform { double x -> diploid.logDensity(x) }

        lrs.@names = std.@names
        
        return lrs
    }

    private logStats(String sample, String chr, Matrix approxStd, List<FastNormal> deletions, FastNormal diploid) {
        log.info "Deletion (mean, stddev) for $sample in $chr at $debugRegionIndex: " + [
            deletions[debugRegionIndex].mean,
            deletions[debugRegionIndex].standardDeviation
        ].join(', ')
            
        double diploidDensity = diploid.logDensity(approxStd[((List<String>)approxStd.getProperty('sample')).indexOf(sample)][debugRegionIndex])
        log.info "Diploid density = " + diploidDensity
        
        double deletionDensity = deletions[debugRegionIndex].logDensity(approxStd[((List<String>)approxStd.getProperty('sample')).indexOf(sample)][debugRegionIndex])
        log.info "Deletion density = " + deletionDensity
        log.info "Likelihood ratio = " + (deletionDensity - diploidDensity)
        if((deletionDensity - diploidDensity) > this.deletionCallThreshold) {
            log.info "Deletion call WILL be emitted for $opts.dr"
        }
        else {
            log.info "Deletion call NOT emitted for $opts.dr"
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
}
