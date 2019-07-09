package gngs

import htsjdk.samtools.util.Interval
import org.junit.Test

import groovy.util.logging.Log

@Log
class LiftOver {
    
    String from
    String to
    
    GenomeResource chainFile
    
    htsjdk.samtools.liftover.LiftOver lo 

    public Region liftOver(Region region) throws Exception {
        
        if(chainFile == null) {
            String toGenome = to[0].toUpperCase() + to.substring(1)
            String urlPattern = "http://hgdownload.cse.ucsc.edu/goldenPath/${from}/liftOver/${from}To${toGenome}.over.chain.gz"
            chainFile = new ResourceDownloader(urlPattern).download()
        }
        
        if(lo == null)
            lo = new htsjdk.samtools.liftover.LiftOver(chainFile.path)
        
        Interval result = lo.liftOver(region as Interval)
        
        Region resultRegion = new Region(result.contig, result.start, result.end)
        resultRegion.setProps(region.expandoProperties)
        
        return resultRegion
    }    
    
    @Test
    void testCall() {
        LiftOver lo = new LiftOver(from:'hg19', to: 'hg38')
        Region result = lo.liftOver(new Region("chr1:1,273,059-1,274,582"))
        println result
    }
}
