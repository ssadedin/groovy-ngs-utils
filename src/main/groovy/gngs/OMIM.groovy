package gngs
/*
 *  Groovy NGS Utils - Some simple utilites for processing Next Generation Sequencing data.
 *
 *  Copyright (C) 2014 Simon Sadedin, ssadedin<at>gmail.com
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
import groovy.transform.CompileStatic;
import java.util.regex.Pattern

/**
 * Parser for OMIM genemap files
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class OMIM {
    
    
    String gene
    
    List<String> disorders = []
    
    Map fields 
    
    public static List GENEMAP2_FIELDS = [
        "Chromosome",
        "start",
        "End",
        "cytoband",
        "cytoloc",
        "MimNo",
        "Genes",
        "GeneName",
        "Symbol",
        "EntrezID",
        "EnsemblGeneID",
        "Comments",
        "Phenotypes",
        "MouseGene"
    ]
    
    final static Pattern TAB = ~/\t/

                
    // # Chromosome Genomic Position Start  Genomic Position End    Cyto Location   Computed Cyto Location  Mim Number  Gene Symbols    Gene Name   Approved Symbol Entrez Gene ID  Ensembl Gene ID Comments    Phenotypes  Mouse Gene Symbol/ID                
                
    //@CompileStatic
    static Map<String,OMIM> parse(String fileName) {
        
        Map<String,OMIM> result = [:]
        
        int recordCount=0
        int geneCount = 0
        new File(fileName).eachLine { String line ->
            
            if(line.startsWith('#'))
                return
            
            String [] parts = OMIM.TAB.split(line)
            Map<String,String> entry = [OMIM.GENEMAP2_FIELDS,parts].transpose().collectEntries()
            
//            println parts
//            println GENEMAP2_FIELDS
//            println entry
//            println "-" * 100
            
//            System.exit(0)
            
            List disorders = ((entry.Phenotypes?.tokenize(";"))?:[]) as List
            
            
            List genes = entry.Genes?.tokenize(",")?:[]
            genes*.trim()?.each { gene ->
                OMIM omim = new OMIM(gene:gene, fields:entry)
                omim.disorders = disorders
                result[gene] = omim
            }
            ++recordCount
        }
        System.err.println "Loaded $recordCount records associated with ${result.size()} genes"
        
        return result
    }
}
