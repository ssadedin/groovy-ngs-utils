package gngs.tools

import gngs.ProgressCounter
import gngs.RefGenes
import gngs.Region
import gngs.Regions
import gngs.ToolBase
import gngs.Utils
import graxxia.Stats
import graxxia.ThresholdRange
import graxxia.Thresholder
import groovy.transform.CompileStatic
import groovy.util.logging.Log

import java.text.NumberFormat

import org.apache.commons.math3.stat.descriptive.SummaryStatistics

@CompileStatic
class CoveragePosition {
    public CoveragePosition(final String chr, final int pos, final List<Integer> covs) {
        this.stats = Stats.from(covs);
        this.chr = chr
        this.pos = pos
    }
    String chr
    int pos
    Stats stats
}

@CompileStatic
class SystematicGap {
    public SystematicGap(String chr) {
        super();
        this.chr = chr;
    }
    String chr
    SummaryStatistics stats = new SummaryStatistics()
}

/**
 * 
 * @author Simon Sadedin
 */
@Log
class SystematicGaps extends ToolBase {
    
    NumberFormat numberFormat = java.text.NumberFormat.numberInstance;

    {
        numberFormat.maximumFractionDigits = 2
        numberFormat.minimumFractionDigits = 0
    }
    
    String separator = '\t'
    
    static void main(String[] args) {
        cli('SystematicGaps -threshold <threshold> <bgzip>  ...', args) {
            threshold 'The coverage threshold below which samples are considered to have insufficent coverage', args:1, required: true
            refgene 'Path to refgene database or auto to download (hg19)', args:1
            format 'csv or tsv for comma vs tab separators', args:1, required: false
            padding 'Amount to widen interval of gap when searching for overlapping genes', args:1
        }
    }
    
    @CompileStatic
    CoveragePosition parseCov(final String line) {
      final List<String> parts = line.tokenize('\t')
      final List<Integer> covs = (parts[3..-1]*.toInteger());
      return new CoveragePosition(parts[0], parts[1].toInteger(),covs)
    }

    @Override
    public void run() {
        if(opts.arguments().size() == 0) 
            throw new IllegalArgumentException("Please provide one or more bgzipped coverage files")
            
        if(opts.csv) 
            separator = ','
        else
        if(opts.tsv) 
            separator = '\t'

        CoveragePosition currentCp = null
        RefGenes refgene = null
        if(opts.refgene) {
            if(opts.refgene == 'auto') {
                log.info "Downloading / loading refgene database ..."
                refgene = RefGenes.download()
            }
            else {
                log.info "Loading refgene database from $opts.refgene ..."
                refgene = new RefGenes(new File(opts.refgene))
            }
        }
        
        GapAnnotator annotator = new GapAnnotator(refgene)
        
        Thresholder<CoveragePosition> thr =  new Thresholder<CoveragePosition>()
        
        final double threshold = opts.threshold.toDouble()
        
        thr.threshold { CoveragePosition cp ->
            return cp.stats.mean - cp.stats.standardDeviation < threshold
        }
        .initWith { CoveragePosition cp, Object old -> 
            currentCp = cp
            new SystematicGap(cp.chr)
        }
        .updateWith { SystematicGap gap, CoveragePosition cp ->
            gap.stats.addValue(cp.stats.mean)
        }
           
        ProgressCounter c = new ProgressCounter(withRate:true, extra:{ "${currentCp?.chr?:'-'} ${thr.ranges.size()} gap regions found"})
        
        Utils.reader(opts.arguments()[0]) { Reader r ->
            r.lines()
             .map(this.&parseCov)
             .forEach { CoveragePosition cpos ->
                 int pos = cpos.pos
                 c.count()
                 thr.update(pos, cpos)
             }
            c.end()
        }
        
        int widenGapForGeneSearchBy = opts.widen ?: 50
                
        log.info "Found ${thr.ranges.size()} systematic coverage gaps"
//        Regions gapRegions = new Regions(thr.ranges.collect { new Region(it.value.chr, it)})
                
        List<String> headers = null

        System.out.withWriter { w ->
            thr.ranges.each { r ->
                SystematicGap gap = r.value

                final Region gapRegion = new Region(gap.chr, r.from, r.to)
                final Region searchRegion = gapRegion.widen(widenGapForGeneSearchBy)
                String gene = searchForOverlappingGenes(refgene, searchRegion)
                List<Map> annotations = annotator.annotateGapRegion(gapRegion)
                
                for(anno in annotations) {
                    anno.remove('gene')
                    if(headers == null) {
                        headers = ['chr','start','end','size','gene','min_cov','max_cov','mean_cov'] + anno*.key
                        w.write(headers.join(separator))
                        w.write('\n')
                    }

                    w.write((
                        [
                            gap.chr, r.from, r.to, r.size(), gene, 
                            numberFormat.format(gap.stats.min), 
                            numberFormat.format(gap.stats.max), 
                            numberFormat.format(gap.stats.mean)
                        ]+anno*.value.collect { (it == null || it == 'null')  ? '' : it }).join(separator)
                    )
                    w.write('\n')
                }
            }
        }
            
    }

    @CompileStatic
    private String searchForOverlappingGenes(RefGenes refgene, Region region) {
        String gene = null
        if(refgene) {
            List<String> genes = refgene.getGenes(region)
            if(genes.size()>1) {
                // If there are multiple genes choose the first that has coding sequence
                // overlapping the region
                gene = genes.find { refgene.getExons(it).overlaps(region) }
            }
            if(!gene)
                gene = genes ? genes[0] : ''
        }
        return gene
    }
}
