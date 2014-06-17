/*
 *  Groovy NGS Utils - Some simple utilites for processing Next Generation Sequencing data.
 *
 *  Copyright (C) 2013 Simon Sadedin, ssadedin<at>gmail.com
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

import java.util.List;


/**
 * An individual in a family pedigree
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class Subject {
    
    String id
    
    Sex sex
    
    List<Integer> phenoTypes
    
    List<Relationship> relationships = []
    
    boolean isAffected() {
        phenoTypes.any { it > 0 }
    }
    
    @CompileStatic
    boolean isChild() {
        relationships.any { Relationship r -> r.type.isChild() }
    }
    
    String toString() {
        id + '(' + sex?.name() + ')'
    }
}

/**
 * Models a relationship between two individuals in a pedigree
 * 
 * @author simon.sadedin@mcri.edu.au
 */
enum RelationshipType {
    
    FATHER {
        double transmits(String chr) {
            if(chr == "chrX")
                0
            else
            if(chr == "chrY")
                1.0
            else
                0.5
        }        
    },
    MOTHER {
        double transmits(String chr) {
            if(chr == "chrY")
                0.0
            else
                0.5
        }                
    },
    SON {
        double transmits(String chr) { 0.0 }        
    },
    DAUGHTER {
        double transmits(String chr) { 0.0 }        
    },
    BROTHER {
        double transmits(String chr) { 0.0 }        
    },
    SISTER {
        double transmits(String chr) { 0.0 }        
    },
    SIBLING {
        double transmits(String chr) { 0.0 }        
    }
     
    double transmits(String chr) {
        return 0.5
    }
    
    boolean isChild() {
        this in [SON,DAUGHTER,BROTHER,SISTER,SIBLING]
    }
}

class Relationship {
    
    RelationshipType type
    
    String from
    
    String to
    
    String toString() {
        type.name() + " to " + to
    }
}

/**
 * Represents a set of relationships between individuals in a family
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class Pedigree {
    
    /**
     * Unique identifier for the family
     */
    String id
    
    /**
     * List of subject ids belonging to the family
     */
//    List<String> samples = []
    
    List<Subject> individuals = []
    
    List<Integer> phenoTypes = []
    
//    Map<String,Relationship> relationships = [:]
    
    String toString() {
        "$id $samples"
    }
    
    Subject motherOf(String id) {
        individuals.find { it.relationships.any { it.type == RelationshipType.MOTHER && it.to == id }}
    }
    
    Subject fatherOf(String id) {
        individuals.find { it.relationships.any { it.type == RelationshipType.FATHER && it.to == id }}
    }
     
    static Map<String,Pedigree> parse(String pedFileName) {
        
        Map<String,Pedigree> families = [:]
        
        List<Subject> subjects = new TSV(pedFileName,columnNames:['familyId','id', 'paternalId', 'maternalId', 'sex', 'phenotype']).collect { line ->
//            println "FamilyId = $line.familyId id = $line.id"
            if(!families.containsKey(line.familyId))
                families[line.familyId] = new Pedigree(id:line.familyId)
            Pedigree p = families[line.familyId]
            def sex = line.sex ?: "other"
            Subject s = new Subject(id: line.id, sex: Sex.decode(sex), phenoTypes:[line.phenotype])
            p.individuals.add(s)
            p.samples.add(s.id)
            p.phenoTypes.add(line.phenotype)
        }
        
        return families
    }
    
    @Lazy
    List<String> affected = { samples.grep {phenoTypes[samples.indexOf(it)] > 0 } }()
    
    @Lazy
    List<String> unaffected = { samples.grep {phenoTypes[samples.indexOf(it)] == 0 } }()
    
    
    List<String> getSamples() {
        individuals*.id
    }
    
    void setSamples(List<String> samples) {
        this.individuals = samples.collect { new Subject(id:it) }
    }
    
    void toPed(Writer w) {
        /* PED Format definition:
         Family ID
         Individual ID
         Paternal ID
         Maternal ID
         Sex (1=male; 2=female; other=unknown)
         Phenotype
         */
        individuals.each { subject ->
            w.println([id, subject.id, motherOf(subject.id)?.id?:"", fatherOf(subject.id)?.id?:"", subject.id in affected ? 1 : 0  ].join("\t"))
        }
    }
}
