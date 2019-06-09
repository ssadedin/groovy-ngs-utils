/*
 *  Groovy NGS Utils - Some simple utilites for processing Next Generation Sequencing data.
 *
 *  Copyright (C) 2018 Simon Sadedin, ssadedin<at>gmail.com
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

import java.nio.file.Files

import groovy.transform.CompileStatic
import groovy.util.logging.Log

@CompileStatic
@Log
class GenomeResource {
    
    File path
    
    boolean stripChr
}

/**
 * Support for downloading data from URLs like UCSC with built-in 
 * caching and genome version detection.
 * 
 * @author Simon Sadedin
 * @param <T>   The class to the instantiated from the downloaded resource
 */
@Log
class ResourceDownloader {
    
    public ResourceDownloader(String urlPattern) {
        this.urlPattern = urlPattern;
    }
    String urlPattern
    
    GenomeResource download(String genomeVersion="hg19") {
        
        Map<String,String> genomeMap = [
            "GRCh37" : "hg19",
            "GRCh38" : "hg38"
        ]
        
        boolean stripChr = false
        
        // Map to appropriate UCSC genome and strip chr if necessary
        String ucscGenomeVersion = genomeVersion
        if(genomeVersion in genomeMap) {
            ucscGenomeVersion = genomeMap[genomeVersion]
            stripChr = true
        }
        
        String fileName = urlPattern.tokenize('/')[-1]
        
        File outputFile = new File(fileName)
        List<File> tried = [outputFile]
        
        try {
            if(!outputFile.exists()) { // In local directory?
                File homeRefGene = new File(System.properties['user.home']+'/' + '.' + fileName)
                if(homeRefGene.exists()) // In home directory?
                    outputFile = homeRefGene
                else
                    tried << homeRefGene
            }
            
            if(!outputFile.exists()) {
                String ucscUrl = urlPattern.replace('##genomeVersion##',ucscGenomeVersion)
                outputFile.withOutputStream { outputStream ->
                    new URL(ucscUrl).withInputStream { urlStream ->
                        Files.copy(urlStream, outputStream)
                    }
                }
            }
            return new GenomeResource(path: outputFile, stripChr:stripChr)
        }
        catch(Exception e) {
//            log.severe("Unable to load downloadable resource after checking in these locations: " + tried.join(','))
            outputFile.delete()
            throw e
        }
    }
}
