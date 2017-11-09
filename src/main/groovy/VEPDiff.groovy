import groovy.util.logging.Log

/**
 * Finds differences in VEP annotations for identical variants in two VCFs
 * 
 * @author Simon Sadedin
 */
@Log
class VEPDiff {
    
    VCF vcf1
    
    VCF vcf2
    
    VEPDiff(VCF vcf1, VCF vcf2) {
        this.vcf1 = vcf1
        this.vcf2 = vcf2
    }
    
    VEPDiff(String vcf1, String vcf2) {
        this(VCF.parse(vcf1), VCF.parse(vcf2))
    }
    
    /**
     */
    void run(PrintWriter out, def outputHTML) {
        
        log.info "Processing VCF1 (${vcf1.size()} variants)"
        log.info "Processing VCF2 (${vcf2.size()} variants)"
        
        // Find all the variants in common
        List<Map> diffs = findDifferences()
        
        log.info "Identified ${diffs.size()} variants in common with differing VEP annotations"
        
        out.println([
            'chr',
            'pos',
            'ref',
            'alt',
            'consdiff',
            'v1cons',
            'v2cons',
            'featdiff',
            'v1feat',
            'v2feat'
        ].join('\t'))
        diffs.each { diff ->
            out.println(
                [
                    diff.v1.chr,
                    diff.v1.pos,
                    diff.v1.ref + '/' + diff.v1.alt,
                    (diff.cons as List).join(','),
                    diff.v1.vepInfo*.Consequence.join(','),
                    diff.v2.vepInfo*.Consequence.join(','),
                    (diff.feat as List).join(','),
                    diff.v1.vepInfo*.Feature.join(','),
                    diff.v2.vepInfo*.Feature.join(',')
                ].join('\t')
                )
        }
    }
    
    /**
     * Return a list of maps containing identical variants having different VEP consequences or features
     * between the two VCFs
     */
    List<Map> findDifferences() {
        ProgressCounter counter = new ProgressCounter(withRate:true, withTime:true, timeInterval:5000)
        List result = vcf1.collect { v1 ->
            counter.count()
            Variant v2 = vcf2.find(v1)
            if(!v2)
              return false
                      
            Set v1Cons = (v1.vepInfo*.Consequence as Set)
            Set v2Cons = (v2.vepInfo*.Consequence as Set)
            Set consDiff = (v1Cons - v2Cons) + (v2Cons - v1Cons)
                    
            Set v1Feat = (v1.vepInfo*.Feature as Set)
            Set v2Feat = (v2.vepInfo*.Feature as Set)
            Set featDiff = (v1Feat - v2Feat) + (v2Feat - v1Feat)
                    
            if(consDiff.size() || featDiff.size())
                return [
                    v1: v1,
                    v2: v2,
                    cons: consDiff,
                    feat: featDiff
                ]
                        
        }.grep {
            it 
        }        
        counter.end()
        return result
    }
    
    static void main(String [] args) {
        
        Utils.configureSimpleLogging()
        
        Cli cli = new Cli(usage:"VEPDiff <options>")
        cli.with {
            vcf1 'First VCF to compare', args:1, required: true
            vcf2 'Second VCF to compare', args:2, required: true
            html 'HTML output file', args:1
            tsv 'TSV output file ', args:1
        }
        
        def opts = cli.parse(args)
        if(!opts) {
            System.exit(1)
        }
        
        String outputPath = opts.tsv?:"-"
        
        PrintWriter out = opts.tsv ? new PrintWriter(new File(outputPath), true) : new PrintWriter(System.out, true)
        new VEPDiff(opts.vcf1, opts.vcf2).run(out, opts.html)
        out.close()
    }
}
