/**
 * Helper utility that walks the user through turning a sample info source / file
 * into a PED file with family relationships.
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class SamplesToPed {

    public SamplesToPed() {
    }
    
    static void main(String [] args) {
        if(!args)
            System.err.println "Please provide sample info file to process"
            
        Map<String,SampleInfo> infos = SampleInfo.parse_sample_info(args[0])
        
        List<String> samples = infos.keySet() as List
        
        infos.each { String sample, SampleInfo info ->
            println sample.padRight(20) + " : " + info.sex
        }
        
        Map<String,Pedigree> pedigrees = [:]
        
        // Find the longest common substring
        Reader r = System.in.newReader()
        def group = lcs(samples)
        
        
        Map<String, Subject> subjects = [:]
        
        group.members.each { sample ->
            println "Pedigree for $sample [$group.group]: "
            String pedigree = r.readLine()
            if(pedigree.isEmpty())
                pedigree = group.group
                
            println "Sex for sample $sample [${infos[sample].sex.name()}]:"
            String sex = r.readLine()
            if(sex.isEmpty())
                sex = infos[sample].sex.name()
                
            println "Relationship for $sample [m/f/c/b/s]: "
            String rel = r.readLine()
            
            String phenotypeValue
            int phenotype = -1
            while(!["a","u"].contains(phenotypeValue)) {
                println "Phenotype for $sample [a/u]: "
                phenotypeValue = r.readLine()
                phenotype = ["a":1, "u":-1][phenotypeValue]
            }
            
            Subject s = new Subject(id: line.id, sex: Sex.decode(sex), phenoTypes:[line.phenotype])
        }
        
        /*
        Pedigree p = pedigrees[pedigree]
        if(!p) {
            p = new Pedigree(id:pedigree)
            pedigrees[pedigree] = p
        }
        */
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
}
