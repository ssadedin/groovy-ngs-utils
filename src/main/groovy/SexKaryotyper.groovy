import gngs.Cli
import gngs.CoverageStats
import gngs.Regions
import gngs.SAM
import gngs.Utils

/**
 * An algorithm based on simple heuristics for estimating the sex of a sample from
 * sequencing coverage.
 * <p>
 * 
 * @author simon
 */
class SexKaryotyper implements Runnable {
    
    SAM bam = null
    
    Regions regions = null
    
    boolean progress = true
    
    CoverageStats yCoverage     
    
    CoverageStats xCoverage 
    
    CoverageStats autosomeCoverage 
    
    List<String> autosomeChrs = ["chr1","chr22"]
    
    Sex sex = null
    
    public SexKaryotyper(SAM bam, Regions regions) {
        this.bam = bam
        this.regions = regions.reduce()
    }
    
    /**
     * The primary logic is simply that lack of coverage of the Y
     * chromosome indicates a female sample. So mean coverage < 5x indicates female.
     * In a simple world, we would end with that and conclude that otherwise the sample
     * is male. There are two things that could go wrong, however: we could have very low
     * coverage overall and thus a male with < 5x, or we could have very high coverage
     * overall and thus a female manages to get > 5x coverage even though she has no Y
     * chromosome. 
     * 
     * The remaining logic protects against these scenarios - first, we only classify
     * a female based on chrY coverage if we have at least 30x mean coverage overall.
     * If we have less than that, or if there is coverage on the Y chromsome, we then
     * look at the ratio of the X chromosome coverage and the autosomes. For a male,
     * we expect it to be about half, while for a female we expect it to be about
     * 1. The 0.7 is simply picking a mid point inbetween those.
     * 
     * Note: chr1 and chr22 are used to generate robust estimates of the mean
     * coverage over the autosomes and use that as an estimate of the expected
     * coverage on the X chromosome for females.
     * 
     * It's worth noting that in Nextera captures, the chrY coverage appears to be
     * artificially boosted - it comes out about the same as autosome coverage even
     * though there is only theoretically 1 copy. For this reason the absolute value
     * of chrY coverage is not compared to autosome coverage in the logic above. 
     */
    void run() {
        
        bam.progress = progress
        
        xCoverage = bam.coverageStatistics(regions.grep { it.chr == "chrX" } as Regions)
        yCoverage = bam.coverageStatistics(regions.grep { it.chr == "chrY" } as Regions)
        autosomeCoverage = bam.coverageStatistics(regions.grep { autosomeChrs.contains(it.chr) } as Regions)
//        autosomeCoverage = bam.coverageStatistics(regions.grep { it.chr != "chrX" && it.chr != "chrY" } as Regions)
        
        if(yCoverage.mean < 5 && xCoverage.mean > 30) {
            sex = Sex.FEMALE
        }
        else
        if(xCoverage.mean / autosomeCoverage.mean < 0.7) {
            sex = Sex.MALE
        }
        else
            sex = Sex.OTHER
    }
    
    static void main(String[] args) {
        Cli cli = new Cli()
        cli.with {
            bam "Alignment of sample to karyotype", args:1, required:true
            bed "BED file of regions to base karyotype on (eg: the BED file of regions covered by sequencing, exome, etc)", args:1, required:true
        }
        
        def opts = cli.parse(args)
        if(!opts) {
            System.exit(1)
        }
        
        SAM sam = new SAM(opts.bam)
        BED bed = new BED(opts.bed).load()
        
        SexKaryotyper k = new SexKaryotyper(sam, bed)
        Utils.time("Running Karyotyper") {
            k.run()
        }
        
        println "xCoverage stats = " + k.xCoverage
        println "=" * 80
        println "yCoverage stats = " + k.yCoverage
        println "=" * 80
        println "autoCoverage stats = " + k.autosomeCoverage
        println "=" * 80
        println "Karyotyping result: " + k.sex
    }
}
