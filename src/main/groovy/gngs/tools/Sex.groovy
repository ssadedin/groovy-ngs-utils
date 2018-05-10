package gngs.tools

import gngs.*

/**
 * Estimates the sex of a sample from one or more VCF file and prints to the console.
 * <p>
 * Usage:
 * <pre>gngs Sex &lt;vcf file&gt;</pre>
 * 
 * See {@link gngs.VCF#guessSex(int)} for details of how sex is determined.
 * 
 * @author Simon Sadedin
 */
class Sex extends ToolBase {
    
    void run() {
        for(String vcf in opts.arguments().grep { it.endsWith('.vcf') }) {
            gngs.Sex sex = new VCF(vcf).guessSex()
            println sex
        }
        
        for(String bamPath in opts.arguments().grep { it.endsWith('.bam') }) {
            if(!opts.t)
                throw new IllegalArgumentException('Please provide a target region to ascertain sex from BAM files')
                
            SexKaryotyper kt = new SexKaryotyper(new SAM(bamPath))
            kt.run()
            println kt.sex
        }        
    }
    
    static main(args) {
        cli("Sex [-t <target region>] <vcf file | bam file> <vcf file | bam file> ...", args) {
            t 'Target regions to analyse (required for BAM files)', longOpt: 'target', args:1, required: false
        }
    }
}
