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

/**
 * Parser for OMIM genemap files
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class OMIM {
    
    
    String gene
    
    List<String> disorders = []
    
    public static List GENEMAP2_FIELDS = [
                "Numbering_System",
                "Month_Entered",
                "Day",
                "Year",
                "CytoLocation",
                "Gene",
                "Status",
                "Title",
                "Title2",
                "MIMNumber",
                "Method",
                "Comments",
                "Comments2",
                "Disorders",
                "Disorders2",
                "Disorders3",
                "Mouse",
                "Reference"
                ]

    //@CompileStatic
    static Map<String,OMIM> parse(String fileName) {
        
        Map<String,OMIM> result = [:]
        
        int recordCount=0
        int geneCount = 0
        new File(fileName).eachLine { String line ->
            
            Map<String,String> entry = [GENEMAP2_FIELDS,line.split("\\|")].transpose().collectEntries()
            
            List disorders = ((entry.Disorders?.split(";"))?:[]) as List
            entry.Gene.split(",")*.trim().each { gene ->
                OMIM omim = new OMIM(gene:gene)
                omim.disorders = disorders
                result[gene] = omim
                ++geneCount
            }
            ++recordCount
        }
        System.err.println "Loaded $recordCount records associated with $geneCount genes"
        
        return result
    }
}
