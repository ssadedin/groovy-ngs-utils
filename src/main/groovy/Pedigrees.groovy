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
 */

import static RelationshipType.*

class Pedigrees {

    /**
     * Lookup table to find family by pedigree Id
     */
    Map<String,Pedigree> families = [:]
    
    /**
     * Look up table to find pedigree by subject Ids
     */
    Map<String, Pedigree> subjects = [:]
    
    Pedigrees add(Pedigrees other) {
        new Pedigrees(families: families + other.families, subjects: subjects + other.subjects)
    }
    
    /**
     * Find the largest possible set of samples that are unrelated to 
     * the given sample. If null is passed, the maximal set of unrelated
     * samples from the entire set will be computed.
     * 
     * @return
     */
    List<Subject> findMaximalUnrelatedSet(String sampleId=null) {
        
        
        List<Pedigree> unrelatedFamilies = null
        
        if(sampleId != null) {
            // Get all the families that are not those of the sample in question
            Pedigree samplePedigree = families[sampleId]
            if(samplePedigree == null) 
                throw new IllegalArgumentException("Sample " + sampleId + " is not a known sample")
                
            unrelatedFamilies = families*.value.grep { it.id != samplePedigree.id}
        }
        else {
            unrelatedFamilies = families*.value
        }
            
        // Each family is independent: choosing the maximal subset of 
        // members from each family will result in the overall maximal number of samples
        List results = []
        for(Pedigree pedigree in unrelatedFamilies) {
            List familyBestSet = pedigree.findMaximalUnrelatedSet()
            results.addAll(familyBestSet)
        }
        return results
    }
    
    /**
     * Parse a PED file to create a set of pedigrees from it
     * 
     * @param pedFileName
     * @return
     */
    static Pedigrees parse(String pedFileName) {
        new File(pedFileName).withReader { r ->
            parse(r)
        }
    }
    
    final static Set<String> NULL_VALUES = [".","0"]
    
    static Pedigrees parse(Reader r) {
        Map<String,Pedigree> families = [:]
        Map<String,Pedigree> subjectsToFamilies = [:]
        List<Subject> subjects = new TSV(r,columnNames:['familyId','id', 'paternalId', 'maternalId', 'sex', 'phenotype']).collect { line ->
            
            if(!families.containsKey(line.familyId))
                families[line.familyId] = new Pedigree(id:line.familyId)
                
            Pedigree p = families[line.familyId]
            def sex = line.sex ?: "other"
            Subject s = new Subject(id: line.id, sex: Sex.decode(sex), phenoTypes:[line.phenotype])
            
            RelationshipType childType = s.sex == Sex.FEMALE ? DAUGHTER : SON
            if(line.paternalId && !NULL_VALUES.contains(line.paternalId))
                s.relationships.add(new Relationship(type:childType,from:s.id, to: line.paternalId))
                
            if(line.maternalId && !NULL_VALUES.contains(line.maternalId))
                s.relationships.add(new Relationship(type:childType,from:s.id, to: line.maternalId)) 
            
            p.individuals.add(s)
            p.phenoTypes.add(line.phenotype)
            
            subjectsToFamilies[s.id] = p
            return s
        }
        
        // The relationships at this point are one-way. We would like to have them two way, so we 
        // go through each subject and if they have a parent, fill them in as a child in the parent's
        // record. Additionally, the relationship between siblings is only implicit in a PED file.
        // As a result, we need to fill it in by inferring it from multiple siblings who have the same
        // parents
        for(Subject s in subjects) {
            
            List<String> parentIds = s.relationships.grep { it.type.isChild() }*.to
            for(String parentId in parentIds) {
                if(!subjectsToFamilies[parentId]) {
                    if(!NULL_VALUES.contains(parentId))
                        System.err.println "WARNING: Sample $s.id has a parent not specified in PED file: $parentId" 
                    continue
                }
                Subject parent = subjectsToFamilies[parentId].individuals.find { it.id == parentId }
                RelationshipType parentType = parent.sex == Sex.FEMALE ? MOTHER : FATHER
                parent.relationships.add(new Relationship(from:parentId, to: s.id, type: parentType))
            }
            
            // Now search for any siblings - ie: any other family member with the same parent
            List<Subject> siblings = subjects.grep { sib -> sib.id != s.id && sib.parentIds.any { it in parentIds } }
            for(Subject sibling in siblings) {
                s.relationships.add(new Relationship(from:s.id, to: sibling.id, type: s.sex == Sex.FEMALE ? SISTER : BROTHER ))
            }
        }
        
        Pedigrees result = new Pedigrees(families:families, subjects: subjectsToFamilies)
        result.families.each { id, ped -> ped.validate() }
        return result
    }
    
	Subject motherOf(String id) {
		subjects[id]?.motherOf(id)
	}
    
	Subject fatherOf(String id) {
		subjects[id]?.fatherOf(id)
	}
        
    void removeFamily(String id) {
        for(String subjectId in families[id].individuals) {
            subjects.remove(subjectId);
        }
        families.remove(id);
    }
    
    List<String> getAffected() {
        families.collect { id, p -> p.affected }.sum()
    }
    
    List<String> getUnaffected() {
        families.collect { id, p -> 
            p.unaffected 
        }.sum()
    }
     
    /**
     * A convenience method that creates a set of pedigrees from set of 
     * unrelated singletons. Such a pedigree is not very useful, but it allows
     * code expecting a pedigree to work with unrelated samples.
     */
    static Pedigrees fromSingletons(List<String> sampleIds) {
        Pedigrees pedigrees = new Pedigrees()
        pedigrees.families = sampleIds.collectEntries { sampleId ->
            Subject s = new Subject(id: sampleId, sex:Sex.UNKNOWN, phenoTypes:[0])
            Pedigree p = new Pedigree(id: sampleId, individuals: [s])
            pedigrees.subjects[sampleId] = p
            [sampleId, p] 
        }
        return pedigrees
    }
    
    String toJson() {
        return "{" + this.families.collect { id, p ->
            "'$id'" + ' : ' + p.toJson() + "\n"
        }.join(",") + "}"
    }
    
    void save(String fileName) {
        new File(fileName).withWriter { w -> save(w) }
    }
    
    void save(Writer w) {
        this.families.each { id, ped ->
            ped.toPed(w)
        }
    }
    
	Subject getAt(String id) {
		subjects[id]?.individuals.find { it.id == id }
	}
    
    void renameSubject(String fromId, String toId) {
        Pedigree family = this.subjects[fromId]
        if(!family) 
            throw new IllegalArgumentException("Sample $fromId is not found in this family")
            
        Subject subject = family.individuals.find { it.id == fromId }
        family.renameSubject(fromId, toId)
        this.subjects.remove(fromId)
        this.subjects[toId] = subject

    }
}
