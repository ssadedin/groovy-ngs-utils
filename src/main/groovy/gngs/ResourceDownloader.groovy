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
        
        Map<String,String> grchToUCSC = [
            "GRCh37" : "hg19",
            "GRCh38" : "hg38"
        ]
        
        Map<String,String> ucscToGRCh = grchToUCSC.collectEntries { [ it.value, it.key.replaceAll('GRCh','') ] }
        
        boolean stripChr = false
        
        // Map to appropriate UCSC genome and strip chr if necessary
        String ucscGenomeVersion = genomeVersion
        if(genomeVersion in grchToUCSC) {
            ucscGenomeVersion = grchToUCSC[genomeVersion]
            stripChr = true
        }
        
        String grchVersion
        if(genomeVersion in ucscToGRCh) {
            grchVersion = ucscToGRCh[grchVersion]
        }
        
        String fileName = genomeVersion ? 
            urlPattern.tokenize('/')[-1].replaceAll('\\.(.*)$', '.' + ucscGenomeVersion + '.$1') : 
            urlPattern.tokenize('/')[-1].replaceAll('\\.(.*)$','.$1') 
        
        File outputFile = new File(fileName)
        List<File> tried = [outputFile]
        
        try {
            if(!outputFile.exists()) { // In local directory?
                File homeRefGene = new File(System.properties['user.home']+'/' + fileName)
                if(homeRefGene.exists()) // In home directory?
                    outputFile = homeRefGene
                else
                    tried << homeRefGene
            }
            
            if(!outputFile.exists()) {
                String resourceUrl = urlPattern.replaceAll('##genomeVersion##',ucscGenomeVersion)
                if(grchVersion)
                    resourceUrl = resourceUrl.replaceAll('##grchVersion##',grchVersion)
                new URL(resourceUrl).withInputStream { urlStream ->
                    Files.copy(urlStream, outputFile.toPath())
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
