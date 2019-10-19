package gngs.tools

import gngs.*

/**
 * Simple utility to extract sample id and sex from VCF file, 
 * and then write that information as a PED file. Parents are left
 * undefined.
 * 
 * @author Simon Sadedin
 */
class VCFtoPed extends ToolBase {
    
    List<VCF> vcfs = null
    
    Map<String,VCF> samplesToVCF = [:]

    @Override
    public void run() {
        
        if(!opts.arguments()) {
            parser.usage()
            System.exit(0)
        }
        
        
        vcfs = opts.arguments().collect { 
            new VCF(it)
        }
        
        // Index the VCFs by sample
        for(VCF vcf in vcfs) {
            vcf.samples.each { samplesToVCF[it] = vcf }
        }
        
        // Write out one line per sample, using the sample itself as the family id
        samplesToVCF.each { sample, vcf ->
            writePEDLine(sample)
        }
    }
    
    void writePEDLine(String sample) {
        VCF vcf = samplesToVCF[sample]
        
        gngs.Sex sex = vcf.guessSex(sample) 
        
        // sex == gngs.Sex.MALE ? 1 : 2,
        
        int pedSex = 0
        if(sex == gngs.Sex.MALE) 
            pedSex = 1
        else
        if(sex == gngs.Sex.FEMALE) 
            pedSex = 2
        
        println([
            sample,
            sample,
            0,
            0,
            pedSex,
            2
        ].join('\t'))
    }
    
    static void main(String[] args) {
        String desc = """
        Utility to write out the samples in one or more VCFs as a PED file, with 
        sex informed by the VCF. Requires indexed VCF.
        """
        cli('VCFtoPed',desc, args) {
        }
    }
}
