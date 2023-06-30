package gngs

import static org.junit.Assert.*

import org.junit.Before
import org.junit.BeforeClass
import org.junit.Test

import groovy.util.logging.Log

@Log
class PubmedTest {
    
    private static final boolean testActualLoad = false
    
    static {
        Utils.configureSimpleLogging()
    }

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

    PubMed  pubmed

    @BeforeClass
    static void deleteBook() {
        log.info "Deleting cached book"
        new File('src/test/data/pubmed_cache/29369572').delete()
    }
    
    void initPubMed() {
        pubmed = new PubMed('src/test/data/pubmed_gene_info.test.tsv','src/test/data/gene2pubmed.test.tsv')
        pubmed.load()
        pubmed.cacheDir = 'src/test/data/pubmed_cache' as File
    }
    
    @Test
    public void testLoadBookAgain() {

        // Disable this test by default because we don't want actual Pubmed getting spammed by testing
        if(!testActualLoad)
            return

        
        initPubMed()


        List<PubMedArticle> articles = pubmed.getArticles([29369572])
        
        println "Retrieved ${articles.size()} articles"
        
        assert articles.size() == 1
        
        assert articles[0].id == 29369572
        
        println articles[0].text

        
        assert articles[0].text.contains('TANGO2 deficiency is characterized')
    }

    @Test
    public void testLoadBook() {

        // Disable this test by default because we don't want actual Pubmed getting spammed by testing
        if(!testActualLoad)
            return

        Utils.configureSimpleLogging()

        PubMed pubmed = new PubMed('src/test/data/pubmed_gene_info.test.tsv','src/test/data/gene2pubmed.test.tsv')
        pubmed.load()
        pubmed.cacheDir = 'src/test/data/pubmed_cache' as File

        List<PubMedArticle> articles = pubmed.getArticles([29369572])
        
        println "Retrieved ${articles.size()} articles"
        
        assert articles.size() == 1
        
        assert articles[0].id == 29369572
        
        println articles[0].text

        
        assert articles[0].text.contains('TANGO2 deficiency is characterized')
    }

    
    @Test
    public void testLoadSpecificArticles() {

        // Disable this test by default because we don't want actual Pubmed getting spammed by testing
        if(!testActualLoad)
            return

        initPubMed()

        List<PubMedArticle> articles = pubmed.getArticles([7566098,15342556])
        
        println "Retrieved ${articles.size()} articles"
        
        assert articles.size() == 2
        
        assert articles[0].id == 7566098
    }

    @Test
    public void testLoadGeneArticles() {

        // Disable this test by default because we don't want actual Pubmed getting spammed by testing
        if(!testActualLoad)
            return

        initPubMed()

        List<PubMedArticle> articles = pubmed.getArticles('TANGO2')
        
        println "Retrieved ${articles.size()} articles"
    }
}
