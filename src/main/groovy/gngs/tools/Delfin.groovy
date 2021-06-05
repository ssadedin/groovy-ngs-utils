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
    double deletionCallThreshold = 3
    
    Regions targets
    
    Matrix covs
    
    Map<String,Regions> results = new ConcurrentHashMap()

    static void main(args) {
        cli('-t <target regions bed> -cov <cov1> [-cov <cov2> ...] -o <deletions.tsv>', args) {
            t 'Target region bed file', args:1, required: true, type: File, longOpt: 'targets'
            o 'Output TSV file containing deletions (stdout)', args:1, required: false, type: File
            lr 'Output file to save likelihood ratio estimates in', args:1, required: false, type: File
            s 'Samples to test for CNVs (all)', args:'*', required: false, type:String
            chr 'Process the given chromosome only (false, process all)', args: 1, required: false
            cov 'Coverage depth file, in GATK format, tab separated with a column of read count for each target region', args: '*', type: File, required: true
        }
    }

    @Override
    public void run() {
        
        log.info "Loading $opts.t"
        targets = new BED(opts.t).load(withExtra:true)
        
        log.info "Loaded ${targets.numberOfRanges} regions (${Utils.humanBp(targets.size())} from $opts.t"
        
        covs = Matrix.concat(*opts.covs.collect { loadCoverage(it) })
        
        log.info "Loaded ${covs.size()} coverage data sets"
        
        List<String> chrs = opts.chr ? [opts.chr] : targets.allRanges*.key
        
        log.info "Processing ${chrs.size()} chromosomes separately in parallel"

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
        
        assert subset != null

        log.info "Calculating likelihoods for $chr"

        Matrix lrs = computeLikelihoods(subset)
        
        // Note: important that samples listed in order they are found in cov matrix
        List<String> samples = opts.ss ? covs.sample.findAll { it in opts.ss }: covs.sample 
        
        // The list of samples to call CNVs for
        List<Integer> sampleRows = opts.ss ? covs.findIndexValues { sample in samples } : (0..<covs.rowDimension)
        
        Matrix sampleLRs = lrs[sampleRows]
        
        if(opts.lr) {
            Utils.writer("${opts.lr}.${chr}.tsv").withWriter { w -> 
                sampleLRs.save(w, r:true) 
            }
        }

        log.info "Calculated sample LRS: \n$sampleLRs"
        
        // Identify candidate CNV regions
        def cnvs = lrs[sampleRows]*.findIndexValues { it > deletionCallThreshold }
        
        log.info "CNV lrs are: \n$cnvs"
        
        Regions chrRegions = targets[colIndexes] as Regions
        Regions mergedCNVs =
            [samples, cnvs].transpose().collectMany { String sampleId, List<Integer> regionIndices ->
                log.info "Found ${regionIndices.size()} for $sampleId in $chr"
                int sampleIndex = covs.sample.indexOf(sampleId)
                chrRegions[regionIndices].eachWithIndex { cnv, index ->
                    cnv.sample = sampleId 
                    cnv.lr = lrs[sampleIndex][(int)regionIndices[index]]
                }
            }

        return mergedCNVs
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
    Matrix computeLikelihoods(final Matrix subset) {

        assert subset != null

        // Diploid regions are by definition std dev 1 because the matrix was standardised
        FastNormal diploid = new FastNormal(0, 1)

        // For deletions we must calculate a separate estimated mean from each column because the 
        // standard deviation of each column is different
        List<FastNormal> deletions = subset.columns.collect {  MatrixColumn c ->
            Stats colStats = Stats.from(c)
            new FastNormal(-0.5 * colStats.mean / colStats.standardDeviation, 1) 
        }
        
        Matrix std = subset.standardise()
        
        Matrix lrs = std.transform { double x, int i, int j -> deletions[j].logDensity(x) } - 
                     std.transform { double x -> diploid.logDensity(x) }

        lrs.@names = std.@names
        
        return lrs
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
