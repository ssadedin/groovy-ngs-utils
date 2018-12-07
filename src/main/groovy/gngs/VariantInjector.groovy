package gngs

import groovy.transform.CompileStatic

/**
 * Injects a variant into a FASTQ read based on recognising the sequence corresponding to 
 * the reference sequence for the location of the variant within the read.
 * 
 * @author Simon Sadedin
 */
class VariantInjector {
    
    FASTA reference
    
    Variant variant
    
    String referenceBases
    
    String replacementBases
    
    String referenceBasesRC
    
    String replacementBasesRC
    
    int nReplaced
    
    int total
    
    int padding
    
    double alleleFraction = 0.5
    
    public VariantInjector(FASTA reference, Variant variant) {
        super();
        this.reference = reference;
        this.variant = variant;
        
        this.variant = variant
        Region variantRegion = new Region(variant).widen(10)
        referenceBases = reference.basesAt(variantRegion)
        
        this.referenceBasesRC= FASTA.reverseComplement(referenceBases)
        
        this.computeReplacementBases()
        
    }
    
    /**
     * Computes the corresponding sequence after a reference sequence is mutated by a variant
     * which must be a SNV, insertion or deletion
     */
    @CompileStatic
    void computeReplacementBases() {
        Allele allele = variant.alleles[0]
        if(allele.type == 'SNP') {
            replacementBases = referenceBases[0..9] + variant.alt + referenceBases[(10+variant.alt.size())..-1]
        }
        else
        if(allele.type == 'DEL')
        {
            replacementBases = referenceBases[0..9] + variant.alt + referenceBases[(11+variant.size())..-1]
            
        }
        else
        if(allele.type == 'INS') {
            replacementBases = referenceBases[0..9] + variant.alt + referenceBases[(10+variant.size())..-1]
        }
        
        this.replacementBasesRC= FASTA.reverseComplement(replacementBases)
    }
    
    /**
     * Injects the variant for this VariantInjector into the given reads, and adjusts
     * the base quality string accordingly.
     * 
     * @param r1    read 1
     * @param r2    read 2
     */
    @CompileStatic
    void inject(FASTQRead r1, FASTQRead r2) {
        
        // We try to always inject into both reads or not at all
        if(inject(r1, referenceBases, replacementBases)) {
            inject(r2, referenceBasesRC, replacementBasesRC, false)
        }
        else // not we don't even attempt to inject into R1 RC if we injected into R1
        if(inject(r2, referenceBases, replacementBases)) {
            inject(r1, referenceBasesRC, replacementBasesRC, false)
        }
    }
    
    @CompileStatic
    boolean inject(FASTQRead r, String refBases, String repBases, boolean random=true) {
        int index = r.bases.indexOf(referenceBases)
        if(index<0) 
            return false
        
        ++total
        if(random && Math.random()>alleleFraction)
            return false
            
        ++nReplaced
        r.bases = r.bases.replace(referenceBases,replacementBases)
        if(r.bases.size() < r.quals.size())
            r.quals = r.quals.substring(0, r.bases.size())
        else
        if(r.bases.size() > r.quals.size()) {
            r.bases = r.bases.substring(0, r.quals.size())
        }
        
        r.name = r.name + ":SIM_${variant.chr}_${variant.pos}_${variant.ref}_${variant.ref}"
    }
}
