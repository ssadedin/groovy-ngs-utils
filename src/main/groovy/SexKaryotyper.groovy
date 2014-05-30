

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
