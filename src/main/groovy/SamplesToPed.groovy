/**
 * Helper utility that walks the user through turning a sample info source / file
 * into a PED file with family relationships.
 * 
 * @author simon.sadedin@mcri.edu.au
 */

import static RelationshipType.*

class SamplesToPed {

    public SamplesToPed() {
    }
    
    static void main(String [] args) {
		
        if(!args) {
            System.err.println "Please provide sample info file to process\n"
            System.exit(1)
        }
            
        Map<String,SampleInfo> infos = SampleInfo.parse_sample_info(args[0])
        
        List<String> samples = infos.keySet() as List
        
        Map<String,Pedigree> pedigrees = [:]
	        
		while(samples) {
	        infos.each { String sample, SampleInfo info ->
	            println sample.padRight(20) + " : " + info.sex
	        }
	        
	        // Find the longest common substring
	        Reader r = System.in.newReader()
	        def group = lcs(samples)
	        
	        
	        Map<String, Subject> subjects = [:]
	        
	        Map<String,Relationship> rels = [:]
	        
	        group.members.each { sample ->
	            println "Pedigree for $sample [$group.group]: "
	            String pedigree = r.readLine()
	            if(pedigree.isEmpty())
	                pedigree = group.group
	                
	                
	            String sex = infos[sample].sex.name()
	            if(infos[sample].sex == Sex.UNKNOWN) {
	                println "Sex for sample $sample [${infos[sample].sex.name()}]:"
	                sex = r.readLine()
	                if(sex.isEmpty()) {
	                    sex = infos[sample].sex.name()
	                }
	            }
	                    
	            String rel 
	            String relDefault = defaultRelationship(sample)
	            while(!["m","f","c","s"].contains(rel)) {
	                println "Relationship for $sample [m/f/c/s] ($relDefault): "
	                String relValue = r.readLine()
	                if(relValue)
	                    rel = relValue
	                else
	                    rel = relDefault
	            }
	            
	            String phenotypeValue
	            int phenotype = rel == "c" ? 1 : 0
	            while(!["a","u"].contains(phenotypeValue)) {
	                println "Phenotype for $sample [a/u] ($phenotype): "
	                phenotypeValue = r.readLine()
	                if(phenotypeValue)
	                    phenotype = ["a":2, "u":1][phenotypeValue]
	                else
	                    break
	            }
	            
	            Subject s = new Subject(id: sample, sex: Sex.decode(sex), phenoTypes:[phenotype])
	            subjects[s.id] = s
	            
	            Pedigree p = pedigrees[pedigree]
	            if(!p) {
	                p = new Pedigree(id:pedigree)
	                pedigrees[pedigree] = p
	            }
	            
	            RelationshipType relationship = decodeRelationship(rel,s)
	          
	            s.relationships.add(new Relationship(type:relationship, from:s.id))
	            
	            p.individuals.add(s)
				
				samples.remove(sample)
	        }
		}
	        
        // Now go through each sample and try to identify the relationships based
        // on the information we do know
        new File(args[0]+".ped").withWriter { w ->
            pedigrees.each { String id, Pedigree p ->
                fillRelationships(p)
                p.toPed(w)
            }
        }
		
		println "Wrote " + args[0]+".ped"
    }
    
    /**
     * Fill in all the relationship ids on the assumption that they are 
     * relative to a single affected child.
     * 
     * @param p
     */
    static void fillRelationships(Pedigree p) {
        for(Subject s in p.individuals) {
            if(s.isChild() || s.relationships.find {it.type in [BROTHER,SISTER]}) {
                linkParentRelationship(p, s, MOTHER)
                linkParentRelationship(p, s, FATHER)                    
            }
        }
    }
    
    /**
     * Find the other individual 
     * 
     * @param p
     * @param s
     * @param r
     */
    static void linkParentRelationship(Pedigree p, Subject s, RelationshipType parentType) {
        
        Relationship r = s.relationships.find { it.type.isChild() }
        if(r.to != null) {
            r = new Relationship(from: s.id, type: r.type)
            s.relationships.add(r)
        }
        
        Subject parent = p.individuals.find { it.relationships.any { it.type == parentType }} 
        if(parent) {
            r.to = parent.id
            def parentRel = parent.relationships.find { it.type == parentType && it.to == null }
            if(!parentRel) {
                parentRel = new Relationship(type:parentType, from: parent.id)
                parent.relationships.add(parentRel)
            }
            parentRel.to = s.id
        }
    }
    
    static RelationshipType decodeRelationship(String rel, Subject s) {
        
        String sex = s.sex.name()
        
        RelationshipType relationship
        if(rel=="m") {
            relationship = MOTHER
        }
        else
        if(rel=="f") {
            relationship = FATHER
        }
        else
        if(rel=="c") {
            if(sex == "MALE")
                relationship = SON
            else
            if(sex == "FEMALE")
                relationship = DAUGHTER
            else
                relationship = SON // ?
        }
        else
        if(rel=="s") {
            if(sex == "MALE")
                relationship = BROTHER
            else
            if(sex == "FEMALE")
                relationship = SISTER
            else
                relationship = SIBLING
        }
        return relationship
    }
    
    static Map lcs(List<String> samples) {
        // Longest sample name
        int len = samples*.size().max()
        def largestGroup
        while(len-->0) {
            // Group by substring
            def groups = samples.groupBy { it.size()>=len?it.substring(0,len): it}
            
            // Largest group
             largestGroup = groups.max { it.value.size() }
             if(largestGroup.value.size() > 1)
                 return [ group: largestGroup.key, members: largestGroup.value ]
        }
        return [ group: largestGroup.key, members: largestGroup.value ]
    }
    
   static String defaultRelationship(String sample) {
       if(sample.endsWith("M"))
           return "m"
       if(sample.endsWith("D"))
           return "f"
       if(sample ==~ '.*[0-9]$')
           return "c"
       return "s"
   }
}
