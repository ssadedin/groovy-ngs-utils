package gngs.tools

import org.apache.commons.math3.stat.descriptive.*

import gngs.*
import graxxia.Matrix
import graxxia.Stats
import groovy.transform.CompileStatic
import groovy.util.logging.Log

@Log
class GeneCovStats extends ToolBase {
    
    List<SAM> bams
    
    RefGenes refGene
    
    Regions targetRegions

    @Override
    public void run() {
        
        log.info "Loading RefSeq exons from $opts.r"
        this.refGene = new RefGenes(opts.r)
        
        this.bams = opts.arguments().collect { new SAM(it) }
        
        log.info "Processing ${bams.size()} BAM files (samples = " + this.bams*.samples*.getAt(0).join(',') + ")"
        List<String> geneList = new File(opts.g).readLines()*.trim()
        
        this.targetRegions = new BED(opts.t, withExtra:true).load(withExtra:true).reduce()
        
        log.info "Loaded ${targetRegions.numberOfRanges} target regions, ${Utils.humanBp(targetRegions.size())}"
        
        log.info "Computing CDS regions for ${geneList.size()} genes ..."
        List<Regions> geneRegions = computeExonicRegions(geneList) 
        
        log.info "Calculating gene coverage statistics for total of ${Utils.humanBp(geneRegions*.size().sum())} ..."
        List<Map> allGeneMasterStats = geneRegions.collect { calculateGeneMasterStats(it) }
        
        log.info "Done"
        
        Matrix m = Matrix.fromListMap(allGeneMasterStats)
        System.out.withWriter { w ->
            m.save(w)
        }
//        println m.toMarkdown().toString()
    }
    
    @CompileStatic
    List<Regions> computeExonicRegions(List<String> genes) {
       return (List<Regions>)genes.collect { symbol ->
           this.refGene.getExons(symbol).each { it.setProperty('gene',symbol) }
       }
    }
    
    public Map<String, Object> calculateGeneMasterStats(Regions geneRegions) {
        
        List<CoverageStats> stats = calcGeneStats(geneRegions)
    
        DescriptiveStatistics master_mean_stats = new DescriptiveStatistics(stats*.mean as double[])
        DescriptiveStatistics master_median_stats = new DescriptiveStatistics(stats*.median as double[])
        DescriptiveStatistics master_mean_gt20 = new DescriptiveStatistics(stats*.fractionAbove(20)*.multiply(100) as double[])
    
        double perc_bases_covered = targetRegions.intersect(geneRegions).size() / geneRegions.size()
    
        return [
            'Gene' : geneRegions[0].gene,
            '% CDS bases covered' : 100* perc_bases_covered,
            'Exon mean of mean coverage' :  master_mean_stats.mean,
            'Exon mean coverage sd' :  master_mean_stats.standardDeviation,
            'Exon mean of median coverage' :  master_median_stats.mean,
            'Exon mean of % > 20x' : master_mean_gt20.mean,
            'Exon mean  % sd' : master_mean_gt20.standardDeviation
        ]
    }
    
    List<CoverageStats> calcGeneStats(Regions gene_regions) {
        return groovyx.gpars.GParsPool.withPool(bams.size()) {
            bams.collectParallel { SAM bam ->
                bam.coverageStatistics(gene_regions)
            }
        }
    }
    
    static void main(String[] args) {
        cli('GeneCovStats -g <gene file> -r <refgene> <bam1> <bam2> ....', args) {
            g 'File containing genes to compute statistics for, one per line', required: true, args: 1
            r 'RefGene database to load gene definitions from (download from UCSC)', required: true, args: 1
            t 'Target regions of analysed / captured regions', required: true, args: 1
        }
    }
}
