package gngs

import static org.junit.Assert.*

import org.junit.Test

import gngs.tools.GapAnnotator

class GapAnnotatorTest {
    
    RefGenes refgenes = new RefGenes('src/test/data/refGene.AR.txt.gz')

    GapAnnotator ga = new GapAnnotator(refgenes)

    @Test
    public void 'test gap intersecting exon is annotated'() {
        
        Region gap = new Region('chrX:66909399-66909460')
        
        List annotations = ga.annotateGapRegion(gap)
        
        println annotations
        
        // There are 4 transcripts, but only 1 intersects an exon sequence, so
        // only that transcript should be written out
        assert annotations.size() == 1
        
        def ann = annotations[0]
        assert ann.id == 'AR'
        assert ann.transcript == 'NM_001348063'
        assert ann.cds_distance == 0 // because it intersects a coding exon
    }
    
    @Test
    void 'gap intersecting intron only is annotated with single transcript with min cds_distance'() {
        
        Region gap = new Region('chrX:66,908,973-66,909,037')
        List annotations = ga.annotateGapRegion(gap)
        
        assert annotations.size() == 1
        
        // check the transcript annotated is really the one with closest CDS distance
        assert annotations[0].transcript == 'NM_001348063'
    }

}
