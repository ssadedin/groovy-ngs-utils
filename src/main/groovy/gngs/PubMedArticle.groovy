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

@Canonical
@CompileStatic
class PubMedArticle {
    String publication
    String title
    String authors
    String authorInfo
    String context
    String text
    String doi
    String pmidInfo
    
    @CompileDynamic
    static PubMedArticle fromXML(GPathResult xml) {
        def article = xml.PubmedArticle.MedlineCitation.Article
        return new PubMedArticle(
            title: article.ArticleTitle.text(),
            text: article.Abstract.AbstractText.text(),
        )
    }
}
