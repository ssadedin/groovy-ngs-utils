/*
 *  Groovy NGS Utils - Some simple utilites for processing Next Generation Sequencing data.
 *
 *  Copyright (C) 2018 Simon Sadedin, ssadedin<at>gmail.com
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
package gngs

import java.util.logging.Level

import gngs.coverage.CoverageCalculatorActor
import gngs.coverage.SampleReadCount
import graxxia.IntegerStats
import groovy.util.logging.Log
import htsjdk.samtools.SAMRecord
import static groovyx.gpars.actor.Actors.actor

/**
 * An algorithm based on simple heuristics for estimating the sex of a sample from
 * sequencing coverage.
 * <p>
 * 
 * @author simon
 */
@Log
class SexKaryotyper implements Runnable {
    
    SAM bam = null
    
    Regions regions = null
    
    boolean progress = true
    
    double yCoverage     
    
    double xCoverage 
    
    double autosomeCoverage 
    
    List<String> autosomeChrs = ["chr1","chr22"]
    
    gngs.Sex sex = null
    
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
        
        String sample = bam.samples[0]
        
        Regions yRegions = (this.regions.grep { it.chr == 'chrY' } as Regions)
        Regions xRegions = (this.regions.grep { it.chr == 'chrX' } as Regions)
        Regions autoRegions = (this.regions.grep { it.chr == 'chr1' || it.chr == 'chr22' } as Regions)        
        
        IntegerStats xStats = calcCovStats(xRegions)
        log.info "Coverage from X calc = $xStats.mean"
        
        IntegerStats yStats = calcCovStats(yRegions)
        log.info "Coverage from Y calc = $yStats.mean"
       
        IntegerStats autoStats = calcCovStats(autoRegions)
        log.info "Coverage from autosome calc = $yStats.mean"
  
        yCoverage = yStats.mean
        xCoverage = xStats.mean
        autosomeCoverage = autoStats.mean
        
        log.info "Means for $sample are: " + [chrX: xCoverage, chrY: yCoverage, autosomes: autosomeCoverage]
        
        if(yCoverage < 5 && xCoverage > 30) {
            sex = Sex.FEMALE
        }
        else
        if(xCoverage / autosomeCoverage < 0.75) {
            log.info "chrX / auto = " + (xCoverage / autosomeCoverage)
            sex = Sex.MALE
        }
        else
            sex = Sex.OTHER
    }
    
    IntegerStats calcCovStats(Regions regions)  {
        IntegerStats covStats = new IntegerStats(1000)
        def coverageCounter = new RegulatingActor<SampleReadCount>(10000, 100000) { 
            void process(SampleReadCount c) {
                covStats.addValue(c.reads)
            }
        }
        coverageCounter.start()
        CoverageCalculatorActor.processBAM(bam, regions, coverageCounter,0)
        coverageCounter.sendStop()
        return covStats
    }
    
    
    static void main(String[] args) {
        
        Cli cli = new Cli()
        cli.with {
            bam "Alignment of sample to karyotype", args:1, required:true
            bed "BED file of regions to base karyotype on (eg: the BED file of regions covered by sequencing, exome, etc)", args:1, required:true
            v 'Enable logging'
        }
        
        def opts = cli.parse(args)
        if(!opts) {
            System.exit(1)
        }
        
        SAM sam = new SAM(opts.bam)
        BED bed = new BED(opts.bed).load()
        
        if(opts.v) {
            Utils.configureSimpleLogging()
        }
        else {
            Utils.configureSimpleLogging(Level.SEVERE)
        }
        
        SexKaryotyper k = new SexKaryotyper(sam, bed)
        Utils.time("Running Karyotyper") {
            k.run()
        }
        
        println "xCoverage mean = " + k.xCoverage
        println "=" * 80
        println "yCoverage mean = " + k.yCoverage
        println "=" * 80
        println "autoCoverage mean = " + k.autosomeCoverage
        println "=" * 80
        println "Karyotyping result: " + k.sex
    }
}
