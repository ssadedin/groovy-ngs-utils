package gngs.tools

import gngs.ToolBase
import gngs.VCF
import gngs.PDF
import graxxia.Binner
import groovy.util.logging.Log
import graxxia.Matrix

@Log
class LCSHPlot extends ToolBase {

    @Override
    public void run() {
        
        new File(opts.o).mkdirs()
        
        log.info "Reading ${opts.arguments()[0]} ..."

        VCF vcf = VCF.parse(opts.arguments()[0], log:log)
        List<String> chrs = (1..22).collect { "chr${it}" } + ["chrX"]
        List<String> files = chrs.collect { chr ->
            log.info "Plot $chr ..."
            plot_chr(vcf, chr)
        }
        
        String sample = vcf.samples[0]

        log.info "Create PDF ..."
        new PDF().document("$opts.o/${sample}_lcsh_report.pdf") {
            title("LCSH Report for $sample")
            br()
            files.each {
                img(it)
            }
        }
        
        log.info "Created $opts.o/${sample}_lcsh_report.pdf"
        log.info "Done"
    }
    
    String plot_chr(VCF vcf, String chr) {
        
        def sample = vcf.samples[0]
        
        def vcfChr = vcf.grep { it.chr == chr }
        
        int maxPos = vcfChr.max { it.pos }.pos
        int binCount = 75
    
        def bins = new Binner(binCount,0,maxPos)
        def stats = bins.stats(vcfChr*.pos, vcfChr.collect { it.het ? 1 : 0 })
    
        def m = new Matrix(pos: bins.midPoints, mean:stats*.mean).fillna(0)
    
        def p = new gngs.plot.Plot(title:"Fraction of Heterozygotes in $sample ($chr)", xLabel: "Position in $chr", xBound: [0,maxPos], yBound: [0,1] ) << \
            new gngs.plot.Line(y:m.mean, x:m.pos)
    
        new File("${sample}_het_dist").mkdirs()
        def fileName = "${opts.o}/${sample}_het_dist_${chr}.png"
        p.save(fileName)
        return fileName
    }
    
    static void main(String[] args) {
        cli('LCSHPlot -o <output directory> <vcf file>', 'Plots the fraction of het variants for diagnosing LCSH regions', args) {
            o 'Directory to write output files to', args:1, required: true
        }
    }
}
