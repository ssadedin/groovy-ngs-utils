package gngs.tools

import gngs.*
import groovy.transform.CompileStatic
import groovy.util.logging.Log

@Log
class VariantInfo extends ToolBase {

    @Override
    public void run() {
        
        VCF vcf = new VCF(opts.i)

        String position = opts.arguments()[0]
        String alt
        if(opts.arguments().size() > 1) {
            def altString = opts.arguments()[1]
            alt = altString.tokenize('/')[-1] // allow either just alt or ref/alt syntax
            log.info "Searching for alt allele $alt based on ${altString}"
        }
        
        boolean hasVep = false
        try {
            vcf.getVepColumns()
            hasVep = true 
        }
        catch(IllegalArgumentException e) {
            // noop
        }

        Variant v = searchForVariant(vcf, position, alt)
        if(!v) {
            println "No variant found matching position $position"
            System.exit(1)
        }
        
        println(""); 
        
        println("$v (QUAL=$v.qual) : ${vcf.samples.join(/,/)}"); 

        println(""); 
        
        // vep is a problem; it is way too long
        // lets break it into each consequence
        Map rawInfo = v.info.grep { !(it.key.toLowerCase() in ['ann','csq','vep'])}.collectEntries()

        Utils.table(title: 'INFO', rawInfo.collect { [ Field : it.key, Value: it.value ] }); 
        
        println("");  

        if(hasVep) {
            List<Map> vepInfos = v.vepInfo
            if(vepInfos)
                Utils.table(title: 'VEP', vepInfos)
            
            println("");  
        }

        Utils.table(title: "Genotypes", [vcf.samples,v.genoTypes].transpose().collect { [ sample: it[0] ] + it[1] } ); println("");
    }
    
    @CompileStatic
    Variant searchForVariant(VCF vcf, String target, String alt) {
        def parts = target.tokenize(':')
        String chr = parts[0]
        int pos = parts[1].toInteger()
        try {
            VCFIndex index = vcf.getIndex()

            log.info "Utilising index to find $target"
            Variant result = null
            index.query(chr, pos-1, pos+1) { if(!result && it.pos == pos && (alt==null || (it.alt == alt))) result = it; }
            return result
        }
        catch(FileNotFoundException e) {
            log.info "No index found; performing linear search for $target ..."
            return (Variant)vcf.iterator().find { Variant v -> v.pos == pos && v.chr == chr && (alt==null || (v.alt == alt)) };
        }
    }
    
    static void main(String[] args) {
        cli('Variant -i <vcf> <position>', 'Show information about a variant in a VCF', args) {
            i 'VCF file - if indexed, will be used, otherwise linear search', type: File, args:1, required: true
        }
    }
}
