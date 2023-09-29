/*
 *  Groovy NGS Utils - Groovy support for processing Next Generation Sequencing data.
 *
 *  Copyright (C) 2023 Simon Sadedin, ssadedin<at>gmail.com
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

import graxxia.TSV

import groovy.transform.CompileStatic
import groovy.util.logging.Log
import groovy.xml.slurpersupport.GPathResult

/**
 * Support for querying and parsing metadata about PubMed articles
 * 
 * Include a caching layer to store queried metadata, to better work
 * with the PubMed API's rate limits.
 * 
 * To use this, you will need to download the source databases for 
 * gene mappings and pubmed article index:
 * 
 * - article index: https://ftp.ncbi.nih.gov/gene/DATA/gene2pubmed.gz
 * - gene mappings: https://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz
 * 
 * Since this class only works with a single taxonomy at a time,
 * you may like to filter the gene mappings file to the specific
 * taxonomy you are using, as it will make load times much faster.
 * 
 * Once you have the database, you can construct a Pubmed instance
 * and call the {@link #load()} method to load the files.
 * 
 * <pre>
 * PubMed pubmed = new PubMed("gene_info.gz", "gene2pubmed.gz")
 * pubmed.load()
 * def articles = pubmed.getArticles("TANGO2")
 * 
 * println "There are ${articles.size()} artices about TANGO2
 * </pre>
 * 
 * To avoid requerying the same data repeatedly, it is recommended
 * to set a cache directory:
 * <pre>
 * pubmed.cacheDir = new File("pubmed_cache")
 * </pre>
 * 
 * See {@link gngs.PubMedArticle} for information about the 
 * article data that is returned.
 * 
 * Note that PubMed implements rate limits, so this class will not attempt
 * to make more than 3 queries per second.
 * 
 * @author simon.sadedin
 */
@Log
class PubMed {

    private static final long serialVersionUID = 1L;
    
    int queryBatchSize = 50
    
    /**
     * We default to human but a different ID could be set here
     */
    int taxonomyId = 9606
    
    Map<String, List<Long>> gene2pmid
    
    Map<String, Long> gene_symbol_to_id
    
    Object symbolFileSource

    Object articleFileSource
    
    /**
     * Directory to cache articles
     */
    File cacheDir
    
    /**
     * For tests
     */
    protected PubMed() {
    }

    PubMed(def geneSymbolFile, def gene2pmid) {
        this.symbolFileSource = geneSymbolFile
        this.articleFileSource = gene2pmid
    }
    
    PubMed load() {
        this.gene2pmid =  parseArticleIndex(articleFileSource)
        this.gene_symbol_to_id = parseGeneSymbolFile(symbolFileSource)
        return this
    }
    
    @CompileStatic
    Map<String, List<Long>> parseArticleIndex(def gene2PubMed) {
        (Map<String, List<Long>>)Utils.reader(gene2PubMed) { r ->
            return new TSV(r, columnNames: ['tax_id', 'GeneID', 'PubMed_ID'], readFirstLine:false).grep {
                it['tax_id']  == taxonomyId
            }
            .groupBy {
                it['GeneID']
            }
            .collectEntries {
                [it.key, ((List<Map>)(it.value))*.PubMed_ID]
            }
        }
    }
    
    Map<String, Long> parseGeneSymbolFile(String geneSymbolFile) {
        new TSV(geneSymbolFile)
        .findAll {
            it['#tax_id'] == taxonomyId
        }
        .collectEntries {
            [it.Symbol, it.GeneID]
        }
    }
    
    /**
     * Parse the raw text from a Pubmed efetch result into a Groovy node hierarchy
     * 
     * @param xmlSource raw text from Pubmed efetch API
     * @return
     */
    GPathResult parseAbstractXML(String xmlSource) {
        def parser = new groovy.xml.XmlSlurper()
        parser.setFeature("http://apache.org/xml/features/disallow-doctype-decl", false)
        parser.setFeature("http://apache.org/xml/features/nonvalidating/load-external-dtd", false)
        return parser.parseText(xmlSource)
    }
    
    /**
     * Query the given list of Pubmed ids efficiently from Pubmed
     * 
     * Tries to balance throttling of the API against size of queries by querying only 
     * {@link #queryBatchSize} articles at a time, while throttling the queries to
     * less than 3/second to respect Pubmed query rate limits.
     * 
     * Articles are cached individually, and only the non-cached articles from the
     * given list will actually be queried.
     * 
     * @param pubmedIds the list of (integer) pubmed ids to be queried
     * @return  a list of {@link PubMedArticle} objects representing the given pubmed ids
     */
    List<PubMedArticle> getArticles(final List<Long> pubmedIds) {
        
        // First populate all entries from cache
        Map<Long, PubMedArticle> cacheQueryResults = queryArticlesFromCache(pubmedIds)
        
        // Get the residual ids that were not successfully queried from the cache
        List<Long> residualIds = cacheQueryResults.grep { it.value == null }*.key
        
        // Query the residual ids using the actual PubMed API in batches of queryBatchSize
        int count = 0
        Map<Long, PubMedArticle> queryResults = 
            residualIds.collate(queryBatchSize).collect { List<Long> subsetOfIds ->
                
                if(subsetOfIds.isEmpty())
                    return [:]
            
                log.info("Querying ${subsetOfIds.size()} pubmed ids [${subsetOfIds.join(',')}]")
                String xmlSource
                URL url = new URL("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=${subsetOfIds.join(',')}&rettype=abstract&retmode=xml")

                // Crude throttling to 3 req/sec
                if(count>3) {
                    Thread.sleep(1100)
                    count = 0
                }
                xmlSource = url.text
                
                ++count

                GPathResult xml = parseAbstractXML(xmlSource)
                cachePubmedResult(xml)

                return [subsetOfIds,PubMedArticle.fromXML(xml)]
                            .transpose()
                            .collectEntries()
            }.sum()?:[:]        
        
        return (cacheQueryResults + queryResults)*.value
    }

    /**
     * If the cache directory is set, save the individual articles from this result as 
     * XML files in the cache.
     * 
     * @param xml   top level XML returned from Pubmed efetch API for abstract query
     */
    private void cachePubmedResult(GPathResult xml) {
        
        Closure cacheResult = { GPathResult pma ->
            Long pmId = (pma.MedlineCitation.PMID.text() ?: pma.BookDocument.PMID.text()).toLong()
            String cacheXML = groovy.xml.XmlUtil.serialize(pma)
            File cacheFile = getCacheFile(pmId)
            if(cacheFile) {
                log.info "Caching article $pmId to ${cacheFile.absolutePath}"
                getCacheFile(pmId).text = cacheXML
            }
        }
        xml.PubmedArticle.each(cacheResult)
        xml.PubmedBookArticle.each(cacheResult)
    }
    
    /**
     * Attempt to query the given list of PubMed ids from the cache. If an article cannot be found
     * in the cache, return it in the result with a null entry.
     * 
     * @param pubmedIds the list of (integer) pubmed ids to be queried
     * @return  
     */
    Map<Long, PubMedArticle> queryArticlesFromCache(List<Long> pubmedIds) {
        
        pubmedIds.collectEntries { pubmedId ->
                    File cacheFilePath = getCacheFile(pubmedId)
                    
                    if(!cacheFilePath || !cacheFilePath.exists())
                        return [ pubmedId, null ] 
                        
                    GPathResult xml = parseAbstractXML(cacheFilePath.text)

                    return [ pubmedId, PubMedArticle.fromXML(xml)[0]]
                }        
    }
    
    @CompileStatic
    List<PubMedArticle> getArticles(String symbol) {
        
        if(this.gene2pmid==null)
            throw new IllegalStateException("Articles not loaded: please call the load() method before querying for articles")

        List<Long> gene_pub_ids = gene2pmid[gene_symbol_to_id[symbol]]
        
        return getArticles(gene_pub_ids)
    }
    
    @CompileStatic
    File getCacheFile(Long pubmedId) {
        if(this.cacheDir == null)
            return null

        if(!this.cacheDir.exists())
            this.cacheDir.mkdirs()

        new File(this.cacheDir, String.valueOf(pubmedId))
    }
}
