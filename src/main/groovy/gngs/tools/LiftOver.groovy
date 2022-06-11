package gngs.tools

import gngs.BED
import gngs.Region
import gngs.Regions
import gngs.ToolBase
import groovy.util.logging.Log

/**
 * Simple tool that translates between genome coordinates with automatic
 * download and caching of the chain files.
 * 
 * @author Simon Sadedin
 */
@Log
class LiftOver extends ToolBase {

    @Override
    public void run() {
        gngs.LiftOver lo = new gngs.LiftOver(from: opts.from, to: opts.to)
        
        def formatResult = opts.tabs ? {
           if(it)
               return ([it.chr, it.from, it.to].join("\t")) 
        } : {
           if(it)
               return (it.toString()) 
        }
        
        int count = 0
        if(opts.bed) {
            Regions regions = new BED(opts.bed, withExtra:true).load(withExtra:true)
            regions.each {  region ->
                def result = lo.liftOver(region)
                if(result) {
                    println([result.chr, result.from, result.to, region.extra].join('\t'))
                }
                ++count
            }
        }
        else
        if(opts.arguments().isEmpty()) {
            System.in.eachLine {  line ->
                println formatResult(lo.liftOver(new Region(line.trim())))
                ++count
            }
        }
        else
        for(String regionValue in opts.arguments()) {
            println formatResult(lo.liftOver(new Region(regionValue)))
                ++count
        }
        
        log.info "Lifted over ${count} coordinates."
    }
    
    static void main(String[] args) {
        cli('Liftover -from <hg38|hg19> -to <hg38|g19> <region1> <region2> ...',
            'Lifts over regions provided as arguments or, if no arguments, reads from standard input', args) {
            bed 'Read regions from BED file <arg>', args: 1, required: false, type: File
            tabs 'Set output format to tab separated'
            from 'Genome build of the source coordinates', args: 1, required: true
            to 'Genome build of the target coordinates', args: 2, required: true
        }
    }
}
