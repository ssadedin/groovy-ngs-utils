package gngs

import static org.junit.Assert.*

import org.junit.Test

class PubmedTest {

    @Test
    public void testLoadPubmed() {
        PubMed pubmed = new PubMed()
        
        def genes_to_articles = pubmed.parseArticleIndex('src/test/data/gene2pubmed.test.tsv')
        
        println "There are ${genes_to_articles.size()} genes with articles"
        
        println "There are ${genes_to_articles*.value*.size().sum()} total articles"
        
        assert genes_to_articles.size() == 2
        
        assert 34949099 in genes_to_articles[6331]
    }

    @Test
    public void testLoadGeneMap() {

        PubMed pubmed = new PubMed()
        
        def pubmed_to_symbol = pubmed.parseGeneSymbolFile('src/test/data/pubmed_gene_info.test.tsv')
        
        println pubmed_to_symbol
        
        assert pubmed_to_symbol.size() == 2
        
        assert pubmed_to_symbol.SCN5A == 6331
    }
    
    
    @Test
    public void testLoadArticle() {

        PubMed pubmed = new PubMed('src/test/data/pubmed_gene_info.test.tsv','src/test/data/gene2pubmed.test.tsv')
        pubmed.load()
        pubmed.cachePath = 'src/test/data/pubmed_cache'

        List<PubMedArticle> articles = pubmed.getArticles('TANGO2')
        
        println "Retrieved ${articles.size()} articles"
    }

}
