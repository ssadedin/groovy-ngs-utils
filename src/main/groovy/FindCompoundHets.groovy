
Cli cli = new Cli()
cli.with {
    vcf 'vcf file to find compount hets in', args:1, required:true
    ped 'ped file describing family relationships - compound hets will only be identified for affected individuals with parents specified', args:1, required:true
    famCount 'maximum number of observations of variant across different families allowed before filtering out (default=3)', args:1
    o 'output file', args:1, required: true
    minQual 'Ignore variants below this quality (default=5.0f)', args:1
    maxMaf 'Maximum MAF for the variants to be considered (default=0.05)', args:1
}

PrintStream log = System.out

opts = cli.parse(args)
if(!opts) {
    System.exit(1)
}

MIN_QUAL = opts.minQual ? opts.minQual.toFloat() : 5.0f;

MAX_MAF = opts.maxMaf ? opts.maxMaf.toFloat() : 0.05f

MAX_FAM_COUNT = opts.famCount ? opts.famCount.toInteger() : 3

Pedigrees pedigrees = Pedigrees.parse(opts.ped)

List<String> unaffected = pedigrees.unaffected
List<String> families = pedigrees.families.keySet() as List
List<String> affected = families.grep { f -> 
    pedigrees.families[f].individuals.grep { Subject s -> s.isChild() && s.isAffected() }
}

MIN_CONSEQUENCE = "coding_sequence_variant"

class Trio {
    String family
    Subject child
    Subject mother
    Subject father
}

class AutoRecTrio {
    Trio trio
    Variant variant
}

List<Trio> trios = affected.collect { s -> 
    new Trio(
        family: pedigrees.subjects[s].id,
        child : pedigrees.subjects[s].individuals.find { it.id == s },
        mother : pedigrees.subjects[s].motherOf(s),
        father : pedigrees.subjects[s].fatherOf(s)
    )
}.grep { it.mother != null && it.father != null }

if(trios.isEmpty()) {
    System.err.println """
        No trios consisting of an affected child and mother and father were found in the VCF file.

        Please ensure the VCF file and PED file contain correct sample information and are consistent with
        each other.
    """.stripIndent()
    System.exit(0)
}

File variantFile = new File(opts.vcf)

//Map<Variant,List<Trio>> autoRecTrios = [:]

List<AutoRecTrio> autoRecTrios = []

// Map by family and then gene to the trios
Map<String, Map<String,List<AutoRecTrio> > > familyAutoRecs = [:]

VCF vcf = VCF.parse(opts.vcf, pedigrees) { Variant v ->

    // Anything with low quality not interesting
    if(v.genoTypes.every { it.GQ < 5.0f}) {
        return false
    }
    
    if(v.getMaxVepMaf()>MAX_MAF)
        return false
        
    if(v.pedigrees.size() > MAX_FAM_COUNT)
        return false
        
    // Check each affected child
    List<Trio> autoRec = trios.grep { Trio trio ->
        (v.sampleDosage(trio.child.id) == 1) && 
            (v.sampleDosage(trio.mother.id) + v.sampleDosage(trio.father.id) == 1)
    }
    
    //log.println "$v has ${autoRec.size()} auto recessive families"
    
    autoRec.each { trio ->
        
        def variantGenes = v.getGenes(MIN_CONSEQUENCE)
        if(!familyAutoRecs.containsKey(trio.family)) {
            familyAutoRecs[trio.family] = [:]
        }
        
        Map genes = familyAutoRecs[trio.family]
        variantGenes.each {  gene ->
//            log.println "Indexing gene $gene for variant $v"
            List<AutoRecTrio> geneVariants = genes[gene]
            if(!geneVariants) {
                genes[gene] = geneVariants = []
            }
            geneVariants.add(new AutoRecTrio(trio: trio, variant:v))
        }
    }
    
    // Only keep variants that satisfy at least the autosomal recessive property
    return !autoRec.isEmpty()
}

// We now have a VCF filtered to contain only variants that have autosomal recessive inheritance 
// in at least one family. We need to find pairs of these that belong to the same family 
// AND are in the same gene AND have opposite inheritance from mother and father
int hetId = 1
List<Variant> toOutput = []
Set<Variant> pendingOutput = new HashSet()
familyAutoRecs.each { String family, Map genes ->
    genes.each { String gene, List<AutoRecTrio> variants ->
        // For each variant in this gene, check if it matches up to any other variant
        for(AutoRecTrio v in variants) {
            for(AutoRecTrio v2 in variants) {
                if(v == v2)
                    break
                if(v.variant.sampleDosage(v.trio.mother.id) == v2.variant.sampleDosage(v.trio.father.id)) {
                    println "Variant $v.variant ${v.variant.vepInfo.grep { it.SYMBOL==gene }*.Consequence.join(',')} and $v2.variant ${v2.variant.vepInfo.grep { it.SYMBOL==gene }*.Consequence.join(',')} are compound het pairs for family $v.trio.family in gene $gene"
                        
                    def updateHet = { Variant vt ->
                        String newChet = "$v.trio.family=$hetId"
                        if(vt.info.CHET)
                            vt.info.CHET = (vt.info.CHET.split(",") + [newChet]).join(",")
                        else
                            vt.info.CHET=newChet
                    }
                    ++hetId
                    
                    if(v.variant.pos == 24460593) {
                        println "Variant $v.variant.pos is paired with $v2.variant.pos"
                    }
                        
                    v.variant.update("Compund Heterozygous ID", updateHet)
                    v2.variant.update("Compund Heterozygous ID", updateHet)
                        
                    if(!pendingOutput.contains(v.variant))
                        toOutput.add(v.variant)
                    if(!pendingOutput.contains(v2.variant))
                        toOutput.add(v2.variant)
                            
                    pendingOutput.add(v.variant)
                    pendingOutput.add(v2.variant)
                }
            }
        }
    }
}
    
new File(opts.o).withOutputStream { o ->
    PrintStream p = new PrintStream(o)
    System.err.println "Writing output to " + opts.o
    vcf.printHeader(p)
    toOutput.each {  p.println it.line }
}    