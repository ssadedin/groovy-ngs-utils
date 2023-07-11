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

import groovy.transform.Canonical
import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.xml.slurpersupport.GPathResult

/**
 * Data structure representing key fields found in PubMed articles and books
 * 
 * @author simon.sadedin
 */
@Canonical
@CompileStatic
class PubMedArticle implements Serializable {
    
    private static final long serialVersionUID = 1L;

    Long id
    String publication
    String title
    String authors
    String authorInfo
    String context
    String text
    String doi
    String pmidInfo
    
    /**
     * Convert a Groovy XML result from PubMed to a PubMedArticle object
     * 
     * Parses both articles or books, representing them using the same
     * structure.
     * 
     * @param xml   abstract XML from the PubMed efetch API
     * @return  list of PubMedArticle objects representing articles and books found in the 
     *          given XML
     */
    @CompileDynamic
    static List<PubMedArticle> fromXML(GPathResult xml) {
        
        def topLevelArticles
        if(xml.MedlineCitation.size()>0) {
            topLevelArticles = [xml]
        }
        else
        if(xml.BookDocument.size()>0) {
            return xml.BookDocument.collect { bookXML ->
                new PubMedArticle(
                    id: bookXML.PMID.text().toLong(),
                    text : bookXML.Abstract.AbstractText.text()
                )
            }
        }
        else {
            topLevelArticles = xml.children()
        }

        return topLevelArticles.collect { topLevelArticle ->
            if(topLevelArticle.MedlineCitation.size()>0) {
                def article = topLevelArticle.MedlineCitation.Article
                return new PubMedArticle(
                    id : topLevelArticle.MedlineCitation.PMID.text().toLong(),
                    title: article.ArticleTitle.text(),
                    text: article.Abstract.AbstractText.text(),
                )
            }
            else {
                def bookXML = topLevelArticle.BookDocument
                return new PubMedArticle(
                    id: bookXML.PMID.text().toLong(),
                    text : bookXML.Abstract.AbstractText.text()
                )
            }
        }
    }
}
