package gngs.tools

import java.awt.event.ItemEvent

import gngs.ToolBase
import gngs.VCF
import gngs.Pedigree
import gngs.Sex
import gngs.Subject
import groovy.transform.CompileStatic
import groovy.transform.ToString
import groovy.util.logging.Log

@ToString(includeNames=true)
@CompileStatic
class PedTrio {
    
    String proband
    String father
    String mother
    
    double denovoRate = 0
}

/**
 * A tool to infer and print out a pedigree file from a trio VCF
 * <p>
 * Relies on the VCF being indexed and having sufficient heterozygous and homozygous variants on the 
 * X chromosome to infer sex.
 * 
 * @author Simon Sadedin
 */
@Log
class TrioPed extends ToolBase {

    @Override
    public void run() {
        
        String vcfPath = opts.arguments()[0]
        
        if(vcfPath == null)
            throw new IllegalArgumentException("Please provide a trio VCF file as an argument")
        
        log.info "Reading $vcfPath"
        VCF vcf = VCF.parse(vcfPath)
        
        List samples = vcf.samples
        if(samples.size() != 3)
            throw new IllegalArgumentException("Please provide a trio VCF file as an argument (wrong number of samples)")
  
        Map<String, Sex> sexes = samples.collectEntries { String sampleId ->
            [sampleId, vcf.guessSex(sampleId)]
        }
        
        List unknownSexSamples = sexes*.key.grep { it == Sex.UNKNOWN }
        if(unknownSexSamples) 
            throw new IllegalStateException("Sex was not able to be inferred for one or more samples: " + unknownSexSamples)
        
        List<List> parentCombos = [samples, samples]
            .combinations()
            .grep { s1, s2 -> s1 != s2 }
            .grep { s1, s2 -> sexes[s1] != sexes[s2] }
            .collect { it.sort { sexes[it] }.reverse() }
            .unique()
        
        log.info "Valid combinations for parents are: " + parentCombos
           
        List<PedTrio> trios = parentCombos.collect { String father, String mother ->
            new PedTrio(father:father, mother: mother, proband: samples.find { !(it in [mother,father]) })
        }
        
        if(trios.every { it.denovoRate > 0.1 }) 
            throw new IllegalStateException("No sample had low rate of de novo variants within trio. Unable to construct valid trio")
        
        // For each trio, find the rate of de novo variants in the child
        log.info "Calculating de novo rates ..."
        trios.each { PedTrio trio ->
            trio.denovoRate = vcf.count { (it.sampleDosage(trio.proband) > 0) && (it.sampleDosage(trio.mother)) == 0 && (it.sampleDosage(trio.father) == 0) } / vcf.size()
            log.info "Trio: $trio"
        }
        
        PedTrio best = trios.min { it.denovoRate }
        log.info "Best trio is: " + best
        
        Subject proband  = new Subject(id:best.proband)
        
        // Finally, write it out as a PED file
        Pedigree ped = new Pedigree(id:opts.family?:best.proband, individuals:[
            proband,
            proband.createFather(best.father),
            proband.createMother(best.mother)
        ])
        
        System.out.withWriter { w -> ped.toPed(w) }
    }
    
    static void main(String[] args) {
        cli('TrioPed <vcf file>', args) {
            family 'Family id (defaults to proband id)', args:1, required: false
        }
    }
}
