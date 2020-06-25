package gngs.tools

import gngs.ProgressCounter
import gngs.ToolBase
import gngs.VCF
import gngs.VCFSiteWalker
import graxxia.Matrix
import groovy.util.logging.Log
import htsjdk.variant.variantcontext.Genotype
import htsjdk.variant.variantcontext.VariantContext

@Log
class SampleCorrelationPlot extends ToolBase {

    @Override
    public void run() {
//        VCF merged = null
//        for(vcfPath in opts.arguments()) {
//            log.info "Reading $vcfPath"
//            VCF vcf = VCF.parse(vcfPath)
//            if(merged == null) {
//                merged = vcf
//            }
//            else {
//                merged = merged.merge(vcf)
//            }
//        }
        
        VariantContext lastVariant = null
        ProgressCounter progress = new ProgressCounter(withRate: true, withTime:true, timeInterval: 10000, extra: {
            "chr: $lastVariant.contig pos: $lastVariant.start (${lastVariant.genotypes[0].sampleName})"
        })
        
        List<String> samples = opts.arguments().collect { new VCF(it).samples }.sum()

        List dosages = []
        
        List snpAlleles = []
        
        String chr = opts.chr?:null

        VCFSiteWalker locusWalker = new VCFSiteWalker(opts.arguments())
        locusWalker.walk { List<VariantContext> variants ->
            
            
            lastVariant = variants[0]
            
            if(chr && lastVariant.getContig() != chr)
                return
                
            if(!variants[0].isSNP())
                return

            int [] variantDosages = new int[samples.size()]
            
            variants.each { VariantContext v ->
                v.genotypes.each { Genotype gt ->

                    if(gt.het) {
                        variantDosages[samples.indexOf(gt.sampleName)] = 1
                    }
                    else
                    if(gt.homVar) {
                        variantDosages[samples.indexOf(gt.sampleName)] = 2
                    }
                }
            }
            
            dosages << variantDosages
            
            snpAlleles << lastVariant.contig + ':' +  lastVariant.start + ' ' + lastVariant.reference +  lastVariant.alternateAlleles[0]
           
            progress.count()
        }
        progress.end()
 
        log.info "Computing SNP matrix"
        
        Matrix snps = new Matrix(dosages)
        
        if(opts.osnps) {
            log.info "Saving SNP matrix to $opts.osnps"
            snps.snp = snpAlleles
            snps.@names = samples
            snps.save(opts.osnps)
        }

        log.info "Computing correlations for ${snps.size()} SNPs"

        Matrix corrs = snps.transpose().rowCorrelations
        corrs.@names = samples
        corrs.Sample = samples
       
        log.info "SNP correlations: \n" + corrs.toMarkdown()
        
        if(opts.ocorr) {
            corrs.save(opts.ocorr)
        }
        
        graxxia.Drawing d = new graxxia.Drawing(opts.o, 120*samples.size(), 100*samples.size(), 0, 0, corrs.columnDimension+2, corrs.columnDimension+1)
        
        int yOffset = 1
        int xOff = 1
        
        d.with {
          color('black')
          fontSize(9)
          corrs.@names.eachWithIndex { sample, i ->
              if(i<corrs.columnDimension-1)
                  text(i+xOff, 4.2, sample)
              
              if(i>0)
                  text(0, i+yOffset+0.5, sample)
          }
          corrs.transform { value, x, y ->
              if(x<y) {
                  color((int)(value*255), 0, 0)
                  d.rect(x+xOff,y+yOffset, x+xOff+1, y+1+yOffset)
              }
              0
          }
          save()
        }
        
    }
    
    static main(args) {
        cli('Plot pairwise SNP correlation between samples', args) {
            o 'Name of image to plot', required:true, args:1
            osnps 'Save SNP matrix in TSV format to this file', args:1
            ocorr 'Save correlation matrix in TSV format to this file', args:1
            chr 'Load SNPs from only this chromosome', args:1
        }
    }

}
