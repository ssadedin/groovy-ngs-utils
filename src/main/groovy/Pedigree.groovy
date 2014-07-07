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

import groovy.json.JsonOutput;
import groovy.transform.CompileStatic;

import java.util.List;

import static RelationshipType.*

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
        phenoTypes.any { it > 1 }
    }
    
    @CompileStatic
    boolean isChild() {
        relationships.any { Relationship r -> r.type.isChild() }
    }
    
    @CompileStatic 
    List<String> getParentIds() {
        relationships.grep { Relationship r -> r.type in [SON,DAUGHTER] }.collect { Relationship r -> r.to }
    }
    
    String toString() {
        id + '(' + sex?.name() + ')'
    }
    
    String toJson() {
        JsonOutput.toJson(
            [
                id : id,
                sex : sex.name(),
                pheno: phenoTypes[0],
                rel: relationships
            ]
        )
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
     * List of subjects belonging to the family
     */
    List<Subject> individuals = []
    
    List<Integer> phenoTypes = []
    
    String toString() {
        "$id $samples"
    }
    
    Subject motherOf(String id) {
        individuals.find { it.relationships.any { it.type == RelationshipType.MOTHER && it.to == id }}
    }
    
    Subject fatherOf(String id) {
        individuals.find { it.relationships.any { it.type == RelationshipType.FATHER && it.to == id }}
    }
     
    
    @Lazy
    List<String> affected = { samples.grep {phenoTypes[samples.indexOf(it)] > 1 } }()
    
    @Lazy
    List<String> unaffected = { samples.grep {phenoTypes[samples.indexOf(it)] <= 1 } }()
    
    
    List<String> getSamples() {
        individuals*.id
    }
    
    List<Subject> findMaximalUnrelatedSet() {
        // We assume that the pedigrees themselves are relatively small. In that case
        // we can just brute force every combination of samples that is unrelated.
        // The number of combinations is 2^n where n is the number of samples in the
        // family.
        final int individualCount = individuals.size()
        final int numCombinations = 1 << individualCount
            
        List best = []
        for(int mask=0; mask<numCombinations; ++mask) {
            // Pull out the samples that are in / out
//            println "Test mask " + String.format("%"+individualCount+"s", Integer.toBinaryString(mask))
            List<Subject> included = []
            
            for(int i=0; i<individualCount; ++i) {
                if(mask & (1<<i)) 
                    included.add(individuals[i])
            }
            
            // No point going further if we would not be better than what we
            // already achieved anyway
            if(included.size() < best.size())
                continue
            
            List<String> includedIds = included*.id
            
            // Are any relatives in there?
//            println "Checking included ids: $includedIds"
            def related = included.grep { s1 -> 
                s1.relationships.any { r -> r.to in includedIds } 
            }
            
            // If they are related, ignore this entire configuration
            if(related) 
                continue
            
            if(included.size()>best.size())
                best = included
        }
        return best
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
            w.println(getPedData(subject).join("\t"))
        }
    }
    
    /**
     * Verify that the samples in the pedigree are internally consistent
     * ie: mother must be female, father male, sample cannot be mother of itself, etc.
     */
    void validate() {
       for(Subject s in individuals) {
           if(s.isChild()) {
               Subject mother = motherOf(s.id)
               if(mother.sex == Sex.MALE) 
                   throw new IllegalStateException("Sample $s.id has mother with sex specified as MALE")
                   
               Subject father = fatherOf(s.id)
               if(father.sex == Sex.FEMALE) 
                   throw new IllegalStateException("Sample $s.id has father with sex specified as FEMALE")
                   
               if(s.id == father.id)
                   throw new IllegalStateException("Sample $s.id has self as father")
                   
               if(s.id == mother.id)
                   throw new IllegalStateException("Sample $s.id has self as mother")
                   
               // More checks - grandparents != self?
               // what else?
           }
       } 
    }
    
    String toJson() {
        "[" + individuals.collect { it.toJson() }.join(",\n") + "]"
    }
    
    List<Object> getPedData(Subject subject) {
       [id, subject.id, fatherOf(subject.id)?.id?:"0",motherOf(subject.id)?.id?:"0",  subject.sex == Sex.MALE ? 1 : 2, subject.id in affected ? 2 : 1  ] 
    }
}
