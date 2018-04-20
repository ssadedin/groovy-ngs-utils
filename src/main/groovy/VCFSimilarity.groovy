
import gngs.Cli
import gngs.ProgressCounter
import gngs.VCF
import gngs.Variant
import graxxia.Matrix
import graxxia.Stats

/**
 * A simple measurement of similarity of different samples by counting the number of 
 * shared vs different SNVs. This measurement can be used as a simplistic way to 
 * detect relatedness among a set of samples that are all sequenced together. For 
 * unrelated samples the number of different SNVs will be approximately similar
 * (while still clustering strongly around population). For related samples there 
 * will be significantly more shared SNVs and fewer distinct SNVs. This tool
 * does not replace a more formal estimation of relatedness, and no hard conclusions 
 * should be drawn from its output. Its main beneifit is that it is very simple and 
 * runs very quickly and easily even on a fairly large set of samples, and thus
 * is suitable for use as a QC measure.
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class VCFSimilarity {
    
    List<String> vcfFiles
    
    List<VCF> vcfs
    
    VCFSimilarity(String... vcfFiles) {
        this.vcfFiles = vcfFiles as List
    }
    
    List<String> samples
    
    Map sampleVariantMetrics = [:]
    
    Matrix sharedCounts 
    
    Matrix onlyCounts 
    
    Matrix zscores
    
    Matrix nonSelfMeans
        
    void compute() {
        
        this.vcfs = vcfFiles.collect { new VCF(it) }
        
        // Find all the samples
        samples = vcfs*.samples.flatten().unique()
        
        def makeMatrix = {
            def m = new Matrix(samples.size(), samples.size())
            m.@names = samples
            m.samples = samples
            return m
        }
        
        sharedCounts = makeMatrix()
        onlyCounts = makeMatrix()
        zscores = makeMatrix()
        
        // Find all the combinations of them all
        List sampleCombinations = GroovyCollections.combinations(samples,samples)
        
        HashMap<String, HashSet> sampleVariants = new HashMap()
        
        for(List samplePair in sampleCombinations) {
            
            String s1 = samplePair[0]
            String s2 = samplePair[1]
            
            println "Analyzing $s1 vs $s2"
            
            int index1 = samples.indexOf(s1)
            int index2 = samples.indexOf(s2)
            
            Set variants1 = sampleVariants[s1]
            if(variants1 == null) {
                VCF vcf1 = vcfs.find { it.samples.contains(s1) }
                variants1 = getVariantSet(vcf1, s1)
                sampleVariants[s1] = variants1
            }
            
            Set variants2 = sampleVariants[s2]
            if(variants2 == null) {
                VCF vcf2 = vcfs.find { it.samples.contains(s2) }
                variants2 = getVariantSet(vcf2, s2)
                sampleVariants[s2] = variants2
            }
             
            // Get superset of all variants for the two samples
            Set union = variants1 + variants2
            
//            println "Variants1 = $variants1"
//            println "Variants2 = $variants1"
//            println "Union = $union"
            
            int count1 = 0
            int count2 = 0
            int both = 0
            
            ProgressCounter c = new ProgressCounter()
            for(String v in union) {
                c.count()
                if((variants1.contains(v)) && (variants2.contains(v))) {
                    ++both
                }
                else 
                if(variants1.contains(v)) {
                    ++count1
                }
                else
                if(variants2.contains(v)) {
                    ++count2
                }
                else
                    assert false : "Variant $v must be in one or the other or both!"
            }
            
            sharedCounts[index1][index2] = both
            sharedCounts[index2][index1] = both
            onlyCounts[index1][index2] = count1
            onlyCounts[index2][index1] = count2
        }
        
        // Compute the means of the non-self, other samples for each sample
        nonSelfMeans = sharedCounts.transform {  c, i, j ->
            List indices = (0..(samples.size()-1)).grep { it != i && it != j }
            Stats.mean(sharedCounts[i][indices])
        }
        
        Matrix sds = sharedCounts.transform {  c, i, j ->
            List indices = (0..(samples.size()-1)).grep { it != i && it != j }
            Stats.from(sharedCounts[i][indices] as double[]).standardDeviation
        }
        
        zscores = (sharedCounts - nonSelfMeans) / sds
        
        zscores.@names = samples
        zscores.samples = samples
    }
   
    Set getVariantSet(VCF vcf, String sample) {
        HashSet variants = new HashSet(30000)
        VCF.parse(vcf.fileName) { Variant v ->
            if(v.sampleDosage(sample)>0 && v.type == "SNP") {
                variants.add("$v.chr:$v.pos:$v.alt".toString())
            }
            return false
        }  
        return variants
    }
    
    public static void main(String [] args) {
        Cli cli = new Cli(usage: 'VCFSimilarity [options] <vcf1> <vcf2> ...')
        cli.with {
            s 'Write matrix of counts of shared SNVs to file', args:1
            o 'Write matrix of counts of distinct SNVs to file', args:1
            z 'Write matrix of z scores to file', args:1
        }
        
        def opts = cli.parse(args)
        if(!opts) {
            cli.usage()
            System.exit(1)
        }
        
        if(opts.arguments().empty) {
            System.err.println "Please provide one or more VCF files to compare\n"
            cli.usage()
            System.exit(1)
        }
        
        VCFSimilarity s = new VCFSimilarity(opts.arguments() as String[])
        
        println "Computing similarity ..."
        s.compute()
        
        println " Shared SNPs ".center(80,"=")
        println s.sharedCounts
        println ""
        
        println " Distinct SNPs ".center(80,"=")
        println s.onlyCounts
        println ""
        
        println " Similarity Scores ".center(80,"=")
        println s.zscores
        println ""
        
        println " Possibly Related Samples ".center(80,"=")
        println ""
        
        List related = []
        s.zscores.transform { value, i, j ->
            if(i == j)
                return 0 
            if(value > 3) {
                related.add([s.samples[i], s.samples[j]].sort().join(" <=> "))
            }
            return  value > 3 ? 1 : 0
        }
        
        related.unique().each {
            println "    " + it
        }
        
        if(opts.s) 
            s.sharedCounts.save(opts.s)
            
        if(opts.o) 
            s.onlyCounts.save(opts.o)
            
        if(opts.z) 
            s.zscores.save(opts.z)
            
        System.out.flush()
    }
}
