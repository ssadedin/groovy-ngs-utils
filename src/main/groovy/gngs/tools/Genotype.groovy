package gngs.tools

import gngs.*
import groovy.util.logging.Log

/**
 * A simple tool that measures the status of predefined SNPs in a 
 * a BAM file, where SNPs are provided as a VCF file.
 * 
 * @author Simon Sadedin
 */
@Log
class Genotype extends ToolBase {

    @Override
    public void run() {
        List<SAM> bams = opts.arguments().collect { new SAM(it) }
        VCF vcf = VCF.parse(opts.vcf)
        
        log.info "Computing ${vcf.size()} * ${bams.size()} genotypes"
        
        List<Map> result = vcf.collect { Variant v ->
            [
                Position: v.chr + ':' + v.pos,
            ] + 
            bams.collectEntries { SAM bam ->
                [ 
                    bam.samples[0],
                    bam.genotype(v.chr, v.pos)
                ]
            }
        }
        
        Utils.table(result)
    }
    
    static void main(String[] args) {
        cli('Genotype -snps <vcf> <bam1> [<bam2>...]', 'Simple SNP genotyping for BAM files', args) {
            vcf 'VCF containing SNPs to test', args:1, required: true
        }
    }
    
}
