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
        
        def ann = annotations.find { it.transcript == 'NM_001348063' }
        
        assert ann != null : 'Expected transcript not annotated: NM_001348063'
        
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
        assert annotations[0].transcript in ['NM_001348063','NM_001348068']
    }

    
    @Test
    void testPanelIndexing() {
        
        ga.indexPanels(["foo","bar"], [
            // First panel (foo)
            [
                [
                    'SubClass1' : 'FOOGENE1',
                    'SubClass2' : 'BARGENE1',
                ]
            ],
            
            // Second panel (bar)
            [
                // this represents a row in the file
                [
                    'SubClass1' : 'FOOGENE2',
                    'SubClass2' : 'BARGENE2'
                ],

                // this represents another row in the file
                [
                    'SubClass1' : 'FOOGENE3',
                    'SubClass2' : 'FOOGENE2' // note: overlaps FOOGENE2 with subclass1
                ]
            ]
        ])
        
        assert 'FOOGENE1' in ga.panelGeneMap
        
        assert 'foo' in ga.panelGeneMap.FOOGENE1

        assert 'BARGENE2' in ga.panelGeneMap
        
        assert ga.panelGeneMap.FOOGENE2.bar.contains('SubClass2')
        assert ga.panelGeneMap.FOOGENE2.bar.contains('SubClass1')
    }
    
    @Test
    void testPanelIndexingFromFile() {
        ga = new GapAnnotator(refgenes, ["src/test/data/test_panel.tsv"])
        
        assert 'CACNA1C' in ga.panelGeneMap
        
        assert 'AKAP9' in ga.panelGeneMap
        
        assert ga.panelGeneMap['CALM1']['test_panel'] 

        assert ga.panelGeneMap['CALM1']['test_panel'].contains('Long QT Syndrome (LQT)')

        assert ga.panelGeneMap['KCNJ2']['test_panel'].size() == 2
        
        assert ga.panelGeneMap['KCNJ2']['test_panel'].contains('Catecholaminergic')
        assert ga.panelGeneMap['KCNJ2']['test_panel'].contains('SomeGenes')

        assert ga.panelGeneMap['CACNA2D1']['test_panel'] == ['Brugada Syndrome']
    }
    
    @Test
    void testPanelAnnotation() {
        ga = new GapAnnotator(refgenes, ["src/test/data/test_panel.tsv", "src/test/data/test_panel2.tsv"])
        Region gap = new Region('chrX:66,908,973-66,909,037')
        List annotations = ga.annotateGapRegion(gap)

        assert annotations.size() == 1

        // check the transcript annotated is really the one with closest CDS distance
        assert annotations[0].transcript in ['NM_001348063','NM_001348068']

        assert annotations[0]['test_panel'] == "no"

        assert annotations[0]["test_panel2"] == "Brugada Syndrome,Catecholaminergic"

        // Case where gene is not in any panel
        ga = new GapAnnotator(refgenes, ["src/test/data/test_panel.tsv"])
        annotations = ga.annotateGapRegion(gap)
        assert annotations[0]['test_panel'] == "no"
    }
}
