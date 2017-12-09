package gngs.tools

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics
import org.apache.commons.math3.stat.descriptive.SummaryStatistics

import gngs.*
import groovy.transform.CompileStatic
import groovy.util.logging.Log

/**
 * Estimates the mean coverage of a BAM file quickly through random sampling.
 * <p>
 * Regions for sampling are selected randomly from a provided BED file. Sampling
 * continues until either the standard deviation of the mean estimate (the mean of the 
 * sampled means) falls below a threshold, or until a maximum number of iterations occurs.
 * 
 * @author Simon Sadedin
 */
@Log
class MeanCoverageEstimator {
    
    SAM bam
    
    Regions regions
    
    Random random 
    
    double sdThreshold = 1.0d
    
    int maxRegions = 1000
    
    int minRegions = 10
    
    int padding = 0
    
    boolean verbose = false

    public MeanCoverageEstimator(SAM sam, Regions regions) {
        this.bam = sam;
        this.regions = regions;
    }
    
    @CompileStatic
    double estimate() {
        
        if(!random)
            random = new Random(0)
        
        DescriptiveStatistics meanStats = new DescriptiveStatistics()
        
        DescriptiveStatistics meanWindow = new DescriptiveStatistics(50)
        
        List<Region> regionList = new ArrayList()
        regionList.addAll(regions)
        
        Set<Integer> usedIndices = new HashSet()
        
        // Randomly sample regions until mean stabilizes
        int regionCount = 0
        while(regionCount++<maxRegions)  {
            
            int regionIndex = (int)(random.nextDouble() * regionList.size())
            if(regionIndex in usedIndices)
                continue
                
            usedIndices << regionIndex
            
            Region region = regionList[regionIndex]
            
            double regionMean = this.bam.meanCoverage(region.chr, region.from-padding, region.to+padding)
            meanStats.addValue(regionMean)
            meanWindow.addValue(meanStats.mean)
            
            if(verbose)
                log.info "$region: $regionMean, N = $meanStats.n, mean = $meanStats.mean median="  + meanStats.getPercentile(50) + " sd = " + meanWindow.standardDeviation
            
            if((regionCount > minRegions) && (meanWindow.standardDeviation<sdThreshold)) {
                if(verbose)
                    log.info "Estimate succeeded"
                return meanWindow.mean
            }
        }
        
        log.info "Exceeded ${maxRegions} iterations"
        return meanWindow.mean
    }
    
    static void main(String [] args) {
        
        Utils.configureSimpleLogging()
        
        Cli cli = new Cli(usage: "MeanCoverageEstimator <args>")
        cli.with {
            bam 'BAM file to estimate coverage for', args:1, required:true
            bed 'Regions to estimate coverage over', args:1, required:true
            padding 'Padding to add to regions', args:1, required:false
            sd 'Threshold for standard deviation for convergence of estimate (accuracy)', args:1, required:false
            seed 'Random seed to use in selecting regions', args:1, required:false
            v 'Verbose mode'
        }
        
        OptionAccessor opts = cli.parse(args)
        if(!opts) 
            System.exit(1)
        
        MeanCoverageEstimator mce = new MeanCoverageEstimator(new SAM(opts.bam), new BED(opts.bed).load())
        if(opts.padding) 
            mce.padding = opts.padding.toInteger()
        
        if(opts.v)
            mce.verbose = true
            
        if(opts.sd) 
            mce.sdThreshold = opts.sd.toDouble()
            
        if(opts.seed)
            mce.random = new Random(opts.seed.toLong())
            
        println mce.estimate()
    }
}
