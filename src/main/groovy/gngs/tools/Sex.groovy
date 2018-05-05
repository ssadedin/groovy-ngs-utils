package gngs.tools

import gngs.Cli
import gngs.ToolBase
import gngs.VCF

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
        for(String vcf in opts.arguments()) {
            gngs.Sex sex = new VCF(vcf).guessSex()
            println sex
        }
    }
    
    static main(args) {
        cli("Sex <vcf file>", args) {
        }
    }
}
