package gngs

import htsjdk.samtools.util.Interval
import groovy.transform.CompileStatic
import groovy.util.logging.Log

@Log
class LiftOver {
    
    String from
    String to
    
    GenomeResource chainFile
    
    htsjdk.samtools.liftover.LiftOver lo 

    public Region liftOver(Region region) throws Exception {
        
        init()
        
        Interval result = lo.liftOver(region as Interval)
        
        if(!result)
            return null
            
        Region resultRegion = new Region(result.contig, result.start, result.end)
        resultRegion.setProps(region.expandoProperties)
        
        return resultRegion
    }
    
    final static Map<String,String> genomeMappings = [
        GRCh37: "hg19",
        GRCh38: "hg38",
        hg19: "hg19",
        hg38: "hg38",
    ]

    @CompileStatic
    private void init() {
        if(chainFile == null) {
            
            String to = genomeMappings.get(to, to)
            String from = genomeMappings.get(from, from)
            
            String toGenome = to[0].toUpperCase() + to.substring(1)
            
            String urlPattern = "https://hgdownload.cse.ucsc.edu/goldenPath/${from}/liftOver/${from}To${toGenome}.over.chain.gz"
            chainFile = new ResourceDownloader(urlPattern).download()
        }

        if(lo == null)
            lo = new htsjdk.samtools.liftover.LiftOver(chainFile.path)
    }    
    
    /**
     * Liftover single coordinate
     */
    @CompileStatic
    public int liftOver(String chr, int pos) throws Exception {
        init()
        Interval i = new Interval(chr, pos, pos+1)
        return lo.liftOver(i).start
    }
    
    
//    @Test
    void testCall() {
        LiftOver lo = new LiftOver(from:'hg19', to: 'hg38')
        Region result = lo.liftOver(new Region("chr1:1,273,059-1,274,582"))
        println result
    }
}
