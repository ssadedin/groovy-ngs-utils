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
package gngs

import groovy.transform.CompileStatic
import groovy.util.logging.Log

/**
 * Injects a variant into a FASTQ read based on recognising the sequence corresponding to 
 * the reference sequence for the location of the variant within the read.
 * 
 * @author Simon Sadedin
 */
@Log
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
    
    boolean debugMode = false
    
    int variantOffset = 0
    
    public VariantInjector(FASTA reference, Variant variant, int variantOffset=0) {
        super();
        this.reference = reference;
        this.variant = variant;
        
        this.variant = variant
        
        
        // variant is wider than that
        if(variantOffset == 0) {
            variantOffset = this.variantOffset = Math.max(10,variant.ref.size()+1) // We only really want 10 bases of context, but we have to include more if the
        }
        else
            this.variantOffset = variantOffset
        
        Region variantRegion = new Region(variant).widen(variantOffset)
        referenceBases = reference.basesAt(variantRegion).toUpperCase()
        
        this.referenceBasesRC= FASTA.reverseComplement(referenceBases)
        
        this.computeReplacementBases()
        
        // log.info Align.global(referenceBases, replacementBases).profile.toString()
    }
    
    /**
     * Computes the corresponding sequence after a reference sequence is mutated by a variant
     * which must be a SNV, insertion or deletion
     */
    @CompileStatic
    void computeReplacementBases() {
        this.replacementBases = variant.applyTo(referenceBases,this.variantOffset+1)
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
        // Note we don't even attempt to inject into R1 RC if we injected into R1
        if(inject(r1, referenceBases, replacementBases)) {
            inject(r2, referenceBasesRC, replacementBasesRC, false)
            renameRead(r1,r2)
        }
        else 
        if(inject(r2, referenceBases, replacementBases)) {
            inject(r1, referenceBasesRC, replacementBasesRC, false)
            renameRead(r1,r2)
        }
        else 
        if(inject(r2, referenceBasesRC, replacementBasesRC)) {
            inject(r1, referenceBases, replacementBases, false)
            renameRead(r1,r2)
        }    
        else 
        if(inject(r1, referenceBasesRC, replacementBasesRC)) {
            inject(r2, referenceBases, replacementBases, false)
            renameRead(r1,r2)
        }        
    }
    
    @CompileStatic
    void renameRead(FASTQRead r1, FASTQRead r2) {
        String suffix =":SIM_${variant.chr}_${variant.pos}_${variant.ref}_${variant.ref}"
        r1.name = r1.name + suffix
        r2.name = r2.name + suffix
    }
  
    
    @CompileStatic
    boolean inject(FASTQRead r, String refBases, String repBases, boolean random=true) {
        int index = r.bases.indexOf(refBases)
        if(index<0) {
            if(debugMode) {
                log.info "No match of $refBases to read $r : $r.bases"
                log.info "Alignment:\n" +Align.global(refBases, r.bases).profile.toString()
            }
            return false
        }
        
        ++total
        if(random && Math.random()>alleleFraction && !debugMode)
            return false
            
        ++nReplaced
        r.bases = r.bases.replace(refBases,repBases)
        if(r.bases.size() < r.quals.size())
            r.quals = r.quals.substring(0, r.bases.size())
        else
        if(r.bases.size() > r.quals.size()) {
            r.bases = r.bases.substring(0, r.quals.size())
        }
       
        if(debugMode)
            log.info "injection to read $r : $r.bases\n" + Align.global(refBases, repBases).profile
            
        return true
    }
}
