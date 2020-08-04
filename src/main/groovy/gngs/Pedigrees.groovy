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
 */

import static gngs.RelationshipType.*

import graxxia.TSV
import static gngs.Sex.*

class Pedigrees {
    
    final static List PED_COLUMNS = ['familyId','id', 'paternalId', 'maternalId', 'sex', 'phenotype']

    /**
     * Lookup table to find family by family Id
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
     * Add the pedigree to this Pedigrees object
     * 
     * @param ped
     */
    void add(Pedigree ped) {
        ped.individuals.each { Subject s ->
            subjects[s.id] = ped
        }
        this.families[ped.id] = ped
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
     * See {@link #parse(Reader,Closure)}
     * 
     * @param pedFileName   path to file
     * @param c             filtering closure
     * 
     * @return
     */
    static Pedigrees parse(String pedFileName, Closure c = null) {
        new File(pedFileName).withReader { r ->
            parse(r,c)
        }
    }
    
    final static Set<String> NULL_VALUES = [".","0"]
    
    /**
     * Parse a PED file from the given reader, with optional filtering via the provided closure
     * 
     * @param r the reader to read the PED file from
     * @param c an optional Closure that will be called for each line (see {@link #PED_COLUMNS})
     *          which will cause samples to be skipped from parsing if it returns false
     * 
     * @return  a Pedigrees object representing the families parsed from the PED file
     */
    static Pedigrees parse(Reader r, Closure c = null) {
        Map<String,Pedigree> families = [:]
        Map<String,Pedigree> subjectsToFamilies = [:]
        Map<String,Subject> validSubjects = [:]
        List<Subject> subjects = new TSV(r,columnNames:PED_COLUMNS, columnTypes:[String,String,String,String,Integer,Integer]).findResults { line ->
            
            if(c != null) {
                if(c(line) == false)
                    return
            }
            
            if(line.id in validSubjects)
                throw new IllegalStateException("Sample $line.id is present multiple times in PED file")
                
            
            if(!families.containsKey(line.familyId))
                families[line.familyId] = new Pedigree(id:line.familyId)
                
            Pedigree p = families[line.familyId]
            def sex = line.sex ?: "other"
            Subject s = new Subject(id: line.id, sex: Sex.decode(sex), phenoTypes:[line.phenotype])
            
            RelationshipType childType = s.sex == Sex.FEMALE ? DAUGHTER : SON
            if(line.paternalId && !NULL_VALUES.contains(line.paternalId)) {
                if(validSubjects[line.paternalId]?.sex == FEMALE)
                    throw new IllegalStateException("Sex for $line.paternalId, father of sample $line.id is declared as female")
                s.relationships.add(new Relationship(type:childType,from:s.id, to: line.paternalId))
            }
              
            if(line.maternalId && !NULL_VALUES.contains(line.maternalId)) {
                if(validSubjects[line.maternalId]?.sex == MALE)
                    throw new IllegalStateException("Sex for $line.maternalId, mother of sample $line.id is declared as male") 
                s.relationships.add(new Relationship(type:childType,from:s.id, to: line.maternalId)) 
            }
            
            p.individuals.add(s)
            p.phenoTypes.add(line.phenotype)
            
            subjectsToFamilies[s.id] = p
            validSubjects[s.id] = s
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
        
//        validate()
        
        return result
    }
    
    void validate() {
        for(Subject s in allSubjects()) {
            Subject m = motherOf(s.id)
            if(m?.sex == MALE)
                throw new IllegalStateException("Mother ${m.id} is specified as male for sample $s.id")
                
            Subject dad = fatherOf(s.id)
            if(dad?.sex == FEMALE)
                throw new IllegalStateException("Father ${dad.id} is specified as female for sample $s.id")
        }
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
    
    List<String> getMales() {
        getSubjectsBySex(Sex.MALE)*.id ?: []
    }
     
    List<String> getFemales() {
        getSubjectsBySex(Sex.FEMALE)*.id ?: []
    }
    
    List<String> getSubjectsBySex(Sex sex) {
        families.collect { id, p -> 
            p.individuals.grep { it.sex == sex }
        }.sum() ?: []
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
    
    void renameFamily(String id, String newId) {
        Pedigree ped = this.families.remove(id)
        if(ped == null)
            throw new IllegalArgumentException("Family $id is not known")
            
        ped.id = newId
        this.families[newId] = ped
    }
    
    List<Subject> getAllSubjects() {
        this.families*.value*.individuals.sum()
    }
    
    /**
     * Get Subject by their id
     * 
     * @param id
     * @return
     */
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
