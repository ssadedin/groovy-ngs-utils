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
 * pubmed.cachePath = "pubmed_cache"
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
class PubMed {
    
    /**
     * We default to human but a different ID could be set here
     */
    int taxonomyId = 9606
    
    Map<String, List<Integer>> gene2pmid
    
    Map<String, Integer> gene_symbol_to_id
    
    Object symbolFileSource

    Object articleFileSource
    
    /**
     * Directory to cache articles
     */
    String cachePath
    
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
    Map<String, List<Integer>> parseArticleIndex(def gene2PubMed) {
        (Map<String, List<Integer>>)Utils.reader(gene2PubMed) { r ->
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
    
    Map<String, Integer> parseGeneSymbolFile(String geneSymbolFile) {
        new TSV(geneSymbolFile)
        .findAll {
            it['#tax_id'] == taxonomyId
        }
        .collectEntries {
            [it.Symbol, it.GeneID]
        }
    }
    
    GPathResult parseAbstractXML(String xmlSource) {
        def parser = new groovy.xml.XmlSlurper()
        parser.setFeature("http://apache.org/xml/features/disallow-doctype-decl", false)
        parser.setFeature("http://apache.org/xml/features/nonvalidating/load-external-dtd", false)
        return parser.parseText(xmlSource)
    }
    
    @CompileStatic
    List<PubMedArticle> getArticles(String symbol) {
        
        if(this.gene2pmid==null)
            throw new IllegalStateException("Articles not loaded: please call the load() method before querying for articles")

        List<Integer> gene_pub_ids = gene2pmid[gene_symbol_to_id[symbol]]
        
        int count = 0
        return gene_pub_ids.collect { pubmedId ->
            
            String xmlSource
            URL url = new URL("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=${pubmedId}&rettype=abstract&retmode=xml")
            File cacheFilePath = new File(this.cachePath, String.valueOf(pubmedId))

            // Crude throttling to 3 req/sec
            Closure doQuery = {
                if(count>3) {
                    Thread.sleep(1100)
                    count = 0
                }
                xmlSource = url.text
            }
            
            ++count

            if(cachePath) {
                if(cacheFilePath.exists()) {
                    xmlSource = cacheFilePath.text
                }
                else {
                    doQuery()
                    cacheFilePath.text = xmlSource
                }
            }
            else {
                doQuery()
            }

            GPathResult xml = parseAbstractXML(xmlSource)
            return PubMedArticle.fromXML(xml)
        }
    }
}
