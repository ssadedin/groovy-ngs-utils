package gngs.tools

import gngs.Cli
import gngs.Utils
import gngs.VCF
import gngs.Variant
import gngs.Utils

class VInfo {
    
    static void main(String [] args) {
        Cli cli = new Cli(usage: "VInfo <position> <vcf>")
        
        OptionAccessor opts = cli.parse(args)
        if(!opts) {
            System.exit(1)
        }
        
        if(opts.arguments().size()<2) {
            System.err.println "\nPlease supply a VCF file to display variant info from\n"
            cli.usage()
            System.err.println ""
            System.exit(0)
        }
        
        int pos = opts.arguments()[0].toInteger()
        
        List<VCF> vcfs = opts.arguments()[1..-1].collect { VCF.parse(it) {
            it.pos == pos
        }}
        
        for(VCF vcf in vcfs) {
            for(Variant v in vcf) {
                println ""
                println "$v :"
                println ""
                summarize(v)
                println ""
            }
        }
    }
    
    static void summarize(Variant v) {
        
        Utils.table(
            [''] + v.header.samples,
            [
                ['Copies'] + v.header.samples.collect { v.sampleDosage(it) },
                ['Depth'] + [v.getAlleleDepths(0), v.getAlleleDepths(1)].transpose()*.join('/'),
            ]
        )
    }

}
