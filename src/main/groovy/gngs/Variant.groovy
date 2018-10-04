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

import groovy.json.JsonOutput;
import groovy.transform.CompileStatic
import java.util.regex.Pattern

/**
 * Support for parsing annotations from SnpEff
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class SnpEffInfo {
    
    /**
     * Valid values of SnpEFF impacts
     */
    static EFFECTS_RANKED = [
        'START_LOST',
        'STOP_GAINED',
        'STOP_LOST',
        'START_GAINED',
        'SPLICE_SITE_ACCEPTOR',
        'SPLICE_SITE_DONOR',
        'FRAME_SHIFT',
        'NON_SYNONYMOUS_CODING',
        'SYNONYMOUS_CODING',
        'EXON',
        'INTRON',
        '3_PRIME',
        '5_PRIME',
        'DOWNSTREAM',
        'UPSTREAM',
        'INTRAGENIC'
    ]
    
    static EFFECT_TO_VEP = [
        'START_LOST' : "transcript_ablation",
        'STOP_GAINED': "stop_gained",
        'STOP_LOST' : "stop_lost",
        'START_GAINED' : "initiator_codon_variant",// ?
        'SPLICE_SITE_ACCEPTOR' : "splice_acceptor_variant",
        'SPLICE_SITE_DONOR' : "splice_donor_variant", 
        'NON_SYNONYMOUS_CODING' : "missense_variant",
        'SYNONYMOUS_CODING' : "synonymous_variant",
        'EXON': "non_coding_exon_variant",
        'INTRON' : "intron_variant",
        '3_PRIME' : "3_prime_UTR_variant",
        '5_PRIME' : "5_prime_UTR_variant", 
        'DOWNSTREAM': "downstream_gene_variant",
        'UPSTREAM' : "upstream_gene_variant",
        'INTERGENIC' : "intergenic_variant"
    ]

    /**
     * The kind of impact. One of the values from {@link SnpEffInfo#EFFECTS_RANKED}
     */
    String type
    
    /**
     * The severity of the impact - HIGH, MODIFIER, MODERATE, LOW
     */
    String impact
    
    /**
     * The gene impacted (if any) by the variant
     */
    String gene
    
    /**
     * The transcript impacted (if any) by the variant
     */
    String transcript
    
    /**
     * The original subclause within the info line from which
     */
    String info
    
    /**
     * This is an alias for the SNPEFF 'type' which is useful when trying to make
     * code use both SNPEFF and VEP annotations
     * 
     * @return the type of SNPEFF effect
     */
    String getConsequence() {
        type
    }
    
    boolean isTruncating() {
        int index = EFFECTS_RANKED.indexOf(type)  
        
        return index >= 0 && index < EFFECTS_RANKED.indexOf('NON_SYNONYMOUS_CODING')
    }
    
    String toString() {
        "type=$type,gene=$gene,rank=$impact"
    }
}

class VEPConsequences {
    static List<String> RANKED_CONSEQUENCES = [
        "transcript_ablation",
        "splice_donor_variant",
        "splice_acceptor_variant",
        "stop_gained",
        "frameshift_variant",
        "stop_lost",
        "initiator_codon_variant",
        "inframe_insertion",
        "inframe_deletion",
        "missense_variant",
        "transcript_amplification",
        "splice_region_variant",
        "incomplete_terminal_codon_variant",
        "stop_retained_variant",
        "coding_sequence_variant",
        "synonymous_variant",
        "mature_miRNA_variant",
        "5_prime_UTR_variant",
        "3_prime_UTR_variant",
        "non_coding_exon_variant",
        "nc_transcript_variant",
        "intron_variant",
        "NMD_transcript_variant",
        "non_coding_transcript_variant",
        "upstream_gene_variant",
        "downstream_gene_variant",
        "TFBS_ablation",
        "TFBS_amplification",
        "TF_binding_site_variant",
        "regulatory_region_variant",
        "regulatory_region_ablation",
        "regulatory_region_amplification",
        "feature_elongation",
        "feature_truncation",
        "intergenic_variant"
    ]
    
    static Map<String,String> VEP_IMPACTS = [
        "transcript_ablation":"HIGH",
        "splice_acceptor_variant":"HIGH",
        "splice_donor_variant":"HIGH",
        "stop_gained":"HIGH",
        "frameshift_variant":"HIGH",
        "stop_lost":"HIGH",
        "start_lost":"HIGH",
        "transcript_amplification":"HIGH",
        "inframe_insertion":"MODERATE",
        "inframe_deletion":"MODERATE",
        "missense_variant":"MODERATE",
        "protein_altering_variant":"MODERATE",
        "splice_region_variant":"LOW",
        "incomplete_terminal_codon_variant":"LOW",
        "stop_retained_variant":"LOW",
        "synonymous_variant":"LOW",
        "coding_sequence_variant":"MODIFIER",
        "mature_miRNA_variant":"MODIFIER",
        "5_prime_UTR_variant":"MODIFIER",
        "3_prime_UTR_variant":"MODIFIER",
        "non_coding_transcript_exon_variant":"MODIFIER",
        "intron_variant":"MODIFIER",
        "NMD_transcript_variant":"MODIFIER",
        "non_coding_transcript_variant":"MODIFIER",
        "upstream_gene_variant":"MODIFIER",
        "downstream_gene_variant":"MODIFIER",
        "TFBS_ablation":"MODIFIER",
        "TFBS_amplification":"MODIFIER",
        "TF_binding_site_variant":"MODIFIER",
        "regulatory_region_ablation":"MODERATE",
        "regulatory_region_amplification":"MODIFIER",
        "feature_elongation":"MODIFIER",
        "regulatory_region_variant":"MODIFIER",
        "feature_truncation":"MODIFIER",
        "intergenic_variant":"MODIFIER"
    ]
    
    /**
     * Returns an integer representing the severity of the consequence of a mutation. 
     * Larger integers represent MORE severe consequences.
     */
    @CompileStatic
    static int severityOf(String cons) {
        
        int index = RANKED_CONSEQUENCES.indexOf(cons)
        if(index<0)
            index = RANKED_CONSEQUENCES.size()
        RANKED_CONSEQUENCES.size() - index
    }
    
    /**
     * A mapping of each VEP consequence to a higher level consequence type
     */
    static Map CONSEQUENCE_TYPES = [
        "upstream_gene_variant":                                        "silent",
        "non_coding_exon_variant":                                      "silent",
        "missense_variant":                                             "missense",
        "downstream_gene_variant":                                      "silent",
        "intron_variant":                                               "silent",
        "splice_region_variant":                                        "splice",
        "synonymous_variant":                                           "silent",
        "nc_transcript_variant":                                        "silent",
        "3_prime_UTR_variant":                                          "silent",
        "5_prime_UTR_variant":                                          "silent",
        "splice_donor_variant":                                         "splice",
        "frameshift_variant":                                           "frameshift",
        "inframe_insertion":                                            "missense",
        "stop_gained":                                                  "nonsense",
        "splice_acceptor_variant":                                      "splice",
        "inframe_deletion":                                             "missense",
        "initiator_codon_variant":                                      "silent",
        "mature_miRNA_variant":                                         "silent",
        "stop_lost":                                                    "nonsense",
        "incomplete_terminal_codon_variant":                            "silent",
        "coding_sequence_variant":                                      "missense",
        "stop_retained_variant":                                        "silent",
        "intergenic_variant":                                            "silent"
       ]
}

/**
 * A specific allele in the context of a Variant in a VCF file
 * <p>
 * Each variant (line) in a VCF file may contain multiple alternate alleles.
 * This class represents one alternate allele in a VCF file.
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class Allele {
    
    int index
    
    /**
     * Actual start position of genomic change for this Allele,
     * relative to reference
     */
    int start
    
   /**
    * Actual end position of genomic change for this Allele
    * relative to reference
    */
    int end
    
    String alt
    
    /**
     * INS, DEL or SNP, INV, DUP, DEL
     */
    String type
    
    int size(String ref) {
        (this.alt.size() == ref.size()) ? 1 : Math.abs(this.alt.size() - ref.size())
    }
    
    String toString() { start != end ? "$start-$end $alt ($type)" : "$start $alt ($type)" }
}

/**
 * Represents the genetic state at a specific locus in a genome
 * ie: captures the information at a single line of a VCF file.
 * Note that this can include multiple alleles that may be present
 * at the site, having different start and end position.
 * <p>
 * The {@link Variant} class includes support for accessing the structured
 * annotations and attributes included in the VCF format and popular
 * annotation tools such as VEP, Annovar and SnpEFF. See {@link #getSnpEffInfo()}
 * and {@link #getVepInfo()} for accessors that return parsed annotation
 * information.
 * <p>
 * <b>Important:</b>The {@link Variant} class represents a line in a VCF file and therefore
 * potentially multiple alleles and corresponding annotations, genotypes, etc. However it is 
 * very common to work with normalised VCFs that guarantee only a single allele per site.
 * To make operations on variants more convenient, some methods have both an allele specific
 * and 'default allele' version. The default allele version will operate on the first allele
 * listed in the alleles. These are convenient to use, but you should always keep in mind they
 * ignore any other alleles and could be unsafe to use on non-normalised VCFs.
 * <p>
 * Although {@link Variant} instances can be constructed manually, typically
 * you will query them from {@link VCF} or {@link VCFIndex} classes, and this
 * is required for some operations because the VCF header is required to 
 * parse and understand some of the fields (for example, sample names). Without
 * a header many operations still work, but you may experience {@link IllegalStateException}
 * exceptions if you call methods that depend on the header.
 * <p>
 * One of the most useful operations is to determine the dosage (ie: heterozygosity, number of copies)
 * of a variant in a particular sample. Eg:
 * <pre>
 * VCF vcf = new VCF.parse("test.vcf")
 * 
 * // Find all the variants that affect sample MYSAMPLE in any way
 * List<Variant> mySampleVariants = vcf.grep {  Variant v -> v.sampleDosage("MYSAMPLE")>0 }
 * 
 * // Find all the homozygous variants (for autosomal + female X chromosomes) 
 * List<Variant> mySampleVariants = vcf.grep {  Variant v -> v.sampleDosage("MYSAMPLE") == 2 }
 * ...
 * </pre>
 * 
 * Some support is implemented for understanding pedigrees via the {@link Pedigree} class.
 * If a pedigree is set on the {@link VCF} class from which the {@link Variant} originates, 
 * then one can query variants to find variants that segregate with families or conditions.
 * <p>
 * A Variant also implements the IRange interface and thus you can use
 * them seamlessly with a {@link RegionSource}. For example, to check
 * if a Variant falls within a BED file works simply:
 * <pre>
 * Variant v = ...
 * BED bed = new BED("test.bed").load()
 * if(v in bed) 
 *  println "Variant v falls in the ranges included in the BED file"
 * </pre>
 * 
 * The Variant class supports working with the full genotype data of each sample, and
 * the genotype fields are parsed for you. For example, to check if every sample is below
 * a particular genotype quality:
 * <pre>
 * if(v.genoTypes.every { it.GQ < 5.0f })
 *     println "Quality is low for every sample!"
 * </pre>    
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class Variant implements IRegion {
    
    ///
    /// Common / shared regexes used in parsing
    ///
    final static Pattern COMMA_SPLIT = ~","
    
    final static Pattern AMPERSAND_SPLIT = ~"&"
    
    final static Pattern  PIPE_SPLIT = ~"\\|"
    
    final static Pattern  COLON_SPLIT = ~":"
    
    final static Pattern PIPE_OR_SLASH_SPLIT = ~'[/|]'
    
    final static Pattern TAB_SPLIT = ~'\t'
 
    /**
     * The chromosome on which the variant falls
     */
    String chr
	
    /**
     * The position at which the reference bases ({@link #ref} attribute) starts.
     * Note that this is <i>not</i> necessarily the start of the genetic change represented
     * by the variant.
     */
    int pos
    
    /**
     * The reference allele at the position. Note that as per VCF spec, this is not necessarily
     * the start of the change indicated by the Variant. For the start of actual changes,
     * you need to query the {@link Allele} objects and look at the {@link Allele#start} and
     * {@link Allele#end} attributes.
     */
    String ref
    
    /**
     * The sequence for the first alternate allele. This is a convenience for the common case where 
     * the first alternate allele is the only one of interest. In general you should be aware that
     * more than one alternate Allele represented by a single Variant.
     */
    String alt
    
    /**
     * List of alternate alleles represented by the variant.
     */
    String [] alts
    
    /**
     * ID column as per VCF spec. Usually this contains the rsID (dbSNP) identifer
     */
    String id
    
    /**
     * One of SNP, INS, DEL or SV indicating the type of change represented by the first 
     * alternate allele. 
     */
    String type
    
    /**
     * The whole INFO field as a raw string. Users will generally call {@link #getInfo()} to get
     * this field in parsed form rather than accessing this field directly.
     */
    String info
    
    /**
     * The filter field as loaded from the VCF file.
     */
    String filter
    
    /**
     * The QUAL field as loaded from the VCF file
     */
    float qual
    
    /**
     * The first character of the first alternate allele in byte form. This is a special case
     * optimization to allow for fast access when interoperating with Picard (which wants to see
     * bytes, not char or Strings).
     */
    byte altByte
    
    /**
     * The VCF from which header information will be extracted when required. This enables the variant to 
     * intelligently parse INFO fields and be aware of sample names. If the header is not provided,
     * many functions still work, but some functions will be disabled.
     */
    VCF header

    Map<String,Object> infos
    
    List<SnpEffInfo> snpEffInfo
    
    /**
     * Cached set of dosages for the first allele.
     * <p>
     * <li> 0 = hom ref
     * <li> 1 = het
     * <li> 2 = hom alt
     * 
     * There is one entry in the list for each sample, in the same order
     * as the samples in the VCF header.
     */
    List<Integer> dosages
    
    Set<Pedigree> pedigrees
    
    List<Map<String,Object>> genoTypes
    
    /**
     * The original line from which this variant was parsed
     */
    String line
    
    List genoTypeFields
	
	
	IntRange getRange() {
        if(isSV()) {
            return pos..pos+this.size()
        }
        else
    		return pos..(pos+alts.max { it.size() }.size())
	}
    
    @CompileStatic
    private boolean parseFields(String line) {
        parseFields(line,true)   
    }
    
    /**
     * Parse the given line from a VCF file.
     * 
     * @param line
     */
    @CompileStatic
    private boolean parseFields(String line, boolean ignoreHomRef) {
        
//        def fields = line.split('[\t ]{1,}')
        List<String>  fields = line.tokenize('\t')
        if(ignoreHomRef && fields[4] == '<NON_REF>')
            return false
        
        chr = fields[0]
        id = fields[2]
        ref = fields[3]
        alts = COMMA_SPLIT.split(fields[4])
        alt = alts[0]
        qual = (fields[5] == ".") ? -1.0f : Float.parseFloat(fields[5])
        filter = fields[6]
        info = fields[7]
        
        pos = Integer.parseInt(fields[1])
        parseGenotypes(fields)
        
        setAlt(alt)
        return true
    }
    
    @CompileStatic
    private void parseGenotypes(List<String> fields) {
        if(fields.size() > 8) {
            
          if(genoTypeFields == null)
              genoTypeFields = fields[8].tokenize(':')
          
          // Split the genoTypes field into separate values and parse them out
          genoTypes = fields[9..-1].collect { String gt -> 
              parseGenoTypeFields(gt)
          } 
        }
    }
    
    // Causes strange typecast error (Char => CharSequence at marked line below
//    @CompileStatic
    Map<String,Object> parseGenoTypeFields(String gt) {
        Set numericGTFields = ["DP","GQ"] as Set
        Set numericListFields = ["AD"] as Set
        [genoTypeFields, gt.tokenize(':')].transpose().collectEntries { Object fieldObj ->
              List field = (List)fieldObj
              String key = (String)field[0]
              String value = String.valueOf(field[1])
              if(key in numericGTFields) { 
                  return [key,convertNumericValue(value)] 
              }
              else
              if(key in numericListFields) {
                  return [key,value.tokenize(",").collect { convertNumericValue(it) }]  // strange typecast exc under compilestatic
              }
              else
                  return field
        } 
    }
    
    @CompileStatic
    Number convertNumericValue(String value, Number defaultValue=0I) {
        if(value == '.')
            return defaultValue
            
        if(value.isInteger())
            return value.toInteger()
        
        return value.toFloat()
    }
    
    /**
     * Convert reference and alternate allele strings into a 
     * mutation type, being one of "SNP","INS","DEL","GAIN"
     * or "LOSS", with the latter two representing CNVs.
     * 
     * @param refSeq
     * @param altSeq
     * @return
     */
    @CompileStatic
    String convertType(String refSeq, String altSeq) {
        String result = "SNP"
        if(altSeq == "DUP" || altSeq == "<DUP>")
            result = "GAIN"
        else
        if(altSeq == "DEL" || altSeq == "<DEL>")
            result = "LOSS"
        else
        if(altSeq == "INV" || altSeq == "<INV>")
            result = "INV"
        else
        if(refSeq.size() < altSeq.size())
            result = 'INS'
        else
        if(refSeq.size() > altSeq.size())
            result = 'DEL'
            
        return result
    }
    
    boolean snpEffDirty = false
    
    /**
     * Update using a default description. 
     * <p>
     * Call this method to prepare the Variant for updating and
     * then make your updates within the Closure that you pass. After 
     * the closure exits, other internal fields that are impacted by your changes
     * will be synchronized with the changes you made.
     * <p>
     * Not all fields are supported! see {@link #update(String, Closure)}
     * <p>
     * It's generally bad form to add INFO fields without annotating what you've done.
     * Prefer to use {@link #update(String, Closure)} and pass in a description
     * over this method.
     * 
     * @param c     Closure within updates can be made.
     */
    void update(Closure c) {
        update("Groovy Variant Processing " + (new Date()).toString(),c)
    }
    
    /**
     * Allows various fields to be updated and then synchronises the 
     * rest of the data with those updated fields
     * <p>
     * Call this method to prepare the Variant for updating and
     * then make your updates within the Closure that you pass. After 
     * the closure exits, other internal fields that are impacted by your changes
     * will be synchronized with the changes you made. 
     * <p>
     * Not all fields are supported! see the ones that are set below.
     * <p>
     * The only update to snpEFF information is to remove individual 
     * annotations.
     */
    void update(String desc, Closure c) {
        this.snpEffDirty = false
        if(line == null)
            line = [chr,pos,id?:".",ref,alts.join(","),".","PASS","ADDED"].join("\t")
            
        c(this)
        
        List<String> fields = TAB_SPLIT.split(line) as List
        fields[0] = chr
        fields[1] = String.valueOf(pos)
        fields[3] = ref
        fields[4] = alt?:"."
        fields[5] = String.format("%2.2f",qual)
        fields[6] = filter
        
        if(snpEffDirty) {
            // The only operation we support on the snpEff meta data is to remove some.
            // Hence we just need to rebuild it by concatenating it
            getInfo().EFF = this.snpEffInfo*.info.join(",")
        }
        
        fields[7] = getInfo().collect {k,v -> v != null?"$k=$v":k}.join(';')
        
        if(this.@genoTypes) {
            rebuildGenotypes(fields)
        }
        
        fields[8] = genoTypeFields.join(':')
        
        line = fields.collect { it?:'' }.join('\t')
        
        // Check if we need to add an INFO header line
        if(this.header != null) {
            getInfo().collect { k,v ->
                if(!header.hasInfo(k)) {
                    header.addInfoHeader(k,desc,v)
                }
           }
        }
        this.@info = null
    }

    /**
     * Update the per-sample genotype info fields by rebuilding them
     * from the genoTypes map.
     * 
     * @param fields    the parsed fields of the VCF line
     */
    @CompileStatic
    private void rebuildGenotypes(List fields) {
        for(int i=0; i<this.header.samples.size(); ++i) {
            fields[9+i] = genoTypeFields.collect { fieldName ->
                def field = genoTypes[i][fieldName]
                field instanceof List ? field.join(',') : field
            }.join(':')
        }
    }
    
    /**
     * VCF requires genotype to be the 1st field, but some tools (R, grrrr) write it in whatever order they 
     * feel like. This function corrects it.
     */
    void fixGenoTypeOrder() {
        def fields = line.split('[\t ]{1,}')
        
        def gtFields = fields[8]
        def gtInfo = fields[9]
        if(!gtFields.startsWith('GT:')) {
            def gtInfos = [gtFields.split(':'),gtInfo.split(':')].transpose().collectEntries { it }
            
            fields[8] = 'GT:' + gtInfos.collect { fieldName, fieldValue -> fieldName }.grep { it != 'GT' }.join(':')
            fields[9] = gtInfos['GT']+ ':' + gtInfos.grep { it.key != 'GT' }.collect { it.value }.join(':')
            line = fields.join('\t')
        }
    }
    
    @CompileStatic
    Map<String,Object> getInfo() {
        if(infos == null) {
            infos =  [:]
            if(info != null) {
                for(String s in info.tokenize(';')) {
                    int i = s.indexOf('=')
                    if(i<0) {
                        infos[s]=Boolean.TRUE
                        continue
                    }
                    infos[s.substring(0,i)] = s.substring(i+1)
                }
            }
        }
        return infos
    }
    
    /**
     * Return a list of SnpEffInfo objects, each describing
     * a separate SnpEff effect caused by the variant
     * @return
     */
    List<SnpEffInfo> getSnpEffInfo() {
        // Since we are handing back the snpeff info, 
        // it is possible the user will modify it
        this.snpEffDirty = true
        
        if(snpEffInfo != null)
            return snpEffInfo
            
        String eff = getInfo()['EFF']
        if(!eff)
            return null
            
        def genes = (eff =~ /[A-Z]*\(([^\)]*)\)/).collect { PIPE_SPLIT.split(it[1]+" ")[5] }
        def txs = (eff =~ /[A-Z]*\(([^\)]*)\)/).collect { PIPE_SPLIT.split(it[1]+" ")[8] }
        def effs = (eff =~ /([A-Z0-9_]*)\(([^\)]*)\)/).collect { it[1] }
        def ranks = (eff =~ /[A-Z]*\(([^\)]*)\)/).collect { PIPE_SPLIT.split(it[1]+" ")[0] }
        def infos = COMMA_SPLIT.split(eff)
        
        snpEffInfo = [genes,effs,ranks,infos,txs].transpose().collect { 
            new SnpEffInfo(gene:it[0], type: it[1], impact:it[2], info: it[3], transcript:it[4]) 
        }
    }
   
    List<Map<String,Object>> getVepInfo() {
        if(this.header == null)
            throw new IllegalStateException("Variant must have a header to query VEP information")
            
        def csqs = getInfo().CSQ 
        def vepFields = null
        if(csqs) {
           vepFields = this.header.getVepColumns("CSQ")
        }
        else
        if((csqs = getInfo().ANN)) {
           vepFields = this.header.getVepColumns("ANN")
        }
        
        if(!csqs)
            return []
//            throw new IllegalStateException("This function requires VEP annotations. Please annotate with VEP in order to proceed")
            
        COMMA_SPLIT.split(csqs).collect { csq -> [vepFields,PIPE_SPLIT.split(csq)].transpose().collectEntries() }
    }
    
    /**
     * Return the most impactful SnpEff effect for this variant
     */
    SnpEffInfo getMaxEffect() {
        def info = getSnpEffInfo()
        return info?.min { int r = SnpEffInfo.EFFECTS_RANKED.indexOf(it.impact); r==null?Integer.MAX_VALUE:r }
    }
    
    /**
     * Return list of dosages (number of copies of allele) for 
     * each sample for the first alternate allele.
     * @return
     */
    List<Integer> getDosages() {
        getDosages(0)
    }
    
    List<Integer> getDosages(int alleleIndex) {
        if(dosages != null && alleleIndex == 0)
            return dosages
        
        // TODO: dosages work differently in de novo callers, 
        // they will just supply the base instead of allele 
        // Here is how it was computed in another script    
        /*
          def alleles = s.split(':')[0].split('/')
          if(alleles[0] in ['A','C','T','G']) 
              alleles[0] == alleles[1] ? 2 : 1
          else { .... }
         */
            
        // Actual format is like this:
        //    1/1 : hom alt 1
        //    1/2 : het alt 1 / alt 2
        //    ./. : not called
        // So to find the dosage for allele 1, we need to split the genotype on slash
        // and then count the number of times the requested allele appears.
        List<Integer> result = genoTypes*.GT.collect{PIPE_OR_SLASH_SPLIT.split(it)}*.count { 
            if(!it.isInteger())
                return 0
                
            if(Integer.parseInt(it) == (alleleIndex+1))
                return 1
                
            return 0
        } 
        
        if(alleleIndex == 0)
            dosages = result
            
        return result
    }

    @CompileStatic 
    int sampleDosage(String sampleName) {
        if(this.header == null)
            throw new IllegalStateException("Variant must have a header to query dosage by sample name")
            
        int sampleIndex = this.header.samples.indexOf(sampleName)
        if(sampleIndex < 0)
            throw new IllegalArgumentException("Sample $sampleName not found in VCF. Known samples are $header.samples")
            
        List<Integer> allDosages = getDosages()
        
        return (int)allDosages[sampleIndex]
    }
    
    @CompileStatic 
    int sampleDosage(String sampleName, int alleleIndex) {
        if(this.header == null)
            throw new IllegalStateException("Variant must have a header to query dosage by sample name")
        int sampleIndex = header.samples.indexOf(sampleName)
        List<Integer> allDosages = getDosages(alleleIndex)
        return (int)allDosages[sampleIndex]
    }    
    
    @CompileStatic 
    /**
     * Return a Map of fields from the genotype column corresponding
     * to the given sample. Common fields include:
     * <li>AD - allele depth (list of numeric values)
     * <li>GQ - genotype qualities for each genotype (string value)
     * <li>GT - list of actual genotypes for sample for each allele (string)
     * <li>DP - total depth (integer)
     * @param sampleName
     */
    Map sampleGenoType(String sampleName) {
        if(this.header == null)
            throw new IllegalStateException("Variant must have a header to query genotypes by sample name")
        int sampleIndex = header.samples.indexOf(sampleName)
        return (Map)genoTypes[sampleIndex]
    }
    
    /**
     * Returns true if this variant segregates in a way compatible with the 
     * phenotype in the given pedigree.
     * That is, if any affected sample (with a phenotype) contains the variant,
     * then all the affected samples must. If no affected samples contain 
     * the variant then this method returns true.
     * <p>
     * Note: this does not check that unaffected samples do NOT contain the 
     * variant.
     * 
     * @param p Pedigree to test 
     * @return
     */
    boolean segregatesWith(Pedigree ped) {
        // Nobody with the variant is affected OR everybody with the variant is affected
        !ped.affected.any {sampleDosage(it)} || ped.affected.every {sampleDosage(it)} 
    }
    
    /**
     * Return a list of all the pedigrees that contain this variant
     */
    Set<Pedigree> getPedigrees() {
        if(this.header == null)
            throw new IllegalStateException("Variant must have a header to query pedigrees")
        if(this.@pedigrees == null) {
          def dsgs = getDosages()
          this.@pedigrees = new HashSet()
          for(int i=0; i<dsgs.size(); ++i) {
              if(dsgs[i] > 0) {
                  def ped = header.findPedigreeBySampleIndex(i)   
                  if(!ped) { 
//                      System.err.println("WARNING: no pedigree information found for sample " + header.samples[i])
//                      System.err.println(Arrays.toString(Thread.currentThread().getStackTrace()))
                  }
                  else
                      this.@pedigrees << ped
              }
          }
        }
        return this.@pedigrees
    }
    
    /**
     * Apply this variant to mutate the given sequence at the specified position
     */
    String applyTo(String sequence, int position) {
        switch(this.type) {
            case "SNP":
                return applySNP(sequence, position)
            case "INS":
                return applyINS(sequence, position)
            case "DEL":
                return applyDEL(sequence, position)
            default:
                throw new IllegalStateException("Unknown type of mutation: $type")
        }
    }
    
    String applySNP(String sequence, int position) {
        StringBuilder result = new StringBuilder()
        result.append(sequence.substring(0,position-1))
        result.append(alt)
        result.append(sequence.substring(position))
        return result.toString()
    }
    
    /**
     * Apply this deletion to the specified sequence at the specified position
     * <p>NOTE: this variant MUST be a deletion
     * @param sequence
     * @param position
     * @return given sequence mutated to include this deletion
     */
    String applyDEL(String sequence, int position) {
        if(this.type != "DEL")
            throw new IllegalStateException("Function can only be called for deletions")
            
        // The actual position of the deletion is the base AFTER the pos field
        sequence.substring(0,position) + sequence.substring(position + ref.size() - alt.size())
    }
    
    /**
     * Apply this deletion to the specified sequence at the specified position
     * <p>NOTE: this variant MUST be a deletion
     * @param sequence
     * @param position
     * @return given sequence mutated to include this deletion
     */
    String applyINS(String sequence, int position) {
        if(this.type != "INS")
            throw new IllegalStateException("Function can only be called for insertions")
            
        // The actual position of the insertion is the base AFTER the pos field
        // but the ALT field includes the prior base too, so go 1 backwards
        sequence.substring(0,position-1) + alt + sequence.substring(position)
    }
    
    /**
     * A debug function, displays the sequence difference and returns the 
     * new sequence
     */
    String displaySequenceDifference(String sequence, int position) {
      String updatedSequence = applyTo(sequence,position)
      if(type == "SNP") {
        [sequence,updatedSequence].each { seq ->
          println "${seq.substring(0,position-1)}[${seq.substring(position-1,position)}]${seq.substring(position)}"
        }
      }
      else 
      if(type == "DEL") {
          println sequence
          println "${sequence.substring(0,position)}${'-'*size()}${sequence.substring(position+size())}"
      }
      else 
      if(type == "INS") {
          println "${sequence.substring(0,position)}${'-'*size()}${sequence.substring(position)}"
          println updatedSequence
      }
      return updatedSequence
    }
    
    @CompileStatic
    boolean isSV() {
        return (alt == "DEL" || alt == "<DEL>" || alt == "DUP" || alt == "<DUP>" || alt == "INV" || alt == "<INV>") 
    }
    
    /**
     * Returns the change in size of the genome caused by this variant's default allele.
     * <p>
     * Note: since point mutations / SNVs don't change the size, they have size zero
     * 
     * @return difference in genome size caused by this variant's default allele
     */
    @CompileStatic
    int size() {
        if(isSV()) {
           Object svLen = this.getInfo().SVLEN
           if(!svLen)
               throw new RuntimeException("VCF file contains structural variants but does not have SVLEN information in INFO field")
           
            return Math.abs(Integer.parseInt(String.valueOf(svLen)))
        }
        else
            return Math.abs(ref.size() - alt.size())
    }

    @CompileStatic
    static Variant parse(String line) {
        parse(line,true)
    }
    
    @CompileStatic
    static Variant parse(String line, boolean ignoreNonRef) {
        Variant parsed = new Variant(line:line)
        if(!parsed.parseFields(line))
            return null
        return parsed
    }
    
    String toString() {
        type in ["GAIN","LOSS"] ? "$chr:$pos-${pos+size()}/$type" : "$chr:$pos $ref/$alt"
    }
    
    /**
     * Return true if the given Annovar formattd position and observation
     * match those of this variant
     * <p>
     * Annovar outputs non-standard formatting that makes it difficult to trace
     * an Annovar variant back to it's VCF source. This method implements
     * the tricky logic to compare an Annovar variant to a VCF equivalent
     * and say if they are the same.
     * 
     * @return if the Annovar variant is not the same variant as this VCF variant,
     *         returns 0. If it is the same, returns the number of the alternate
     *         allele (if the variant has only one alternate, the return value will 
     *         be 1).
     */
    @CompileStatic
    int equalsAnnovar(String chr, int pos, String obs) {
        if(this.chr != chr) 
            return 0
            
        int count = 1
        for(alleleTypePair in this.getAllelesAndTypes()) {
            
            String alleleAlt = alleleTypePair[0]
            String alleleType = alleleTypePair[1]
            
            if(this.pos == pos && alleleAlt == obs)  
                return count
              
            // Simple insertion
            if(this.pos == (pos-this.ref.size()+1) && alleleType=="INS" && alleleAlt.endsWith(obs))
                return count
                
            // More complex case - we actually have to splice the insertion in and see if the
            // result matches the reference sequence
            // See VariantTest.testComplexIndelAnnovar
            if(alleleType == "INS") {
                // Modify the reference sequence by adding the insertion
                int offset = pos-this.pos+1
                if(offset <= this.ref.size()) {
                    String newSeq = this.ref.substring(0,offset) + obs + ((offset<this.ref.size()) ? this.ref.substring(offset) : "")
                    if(newSeq==alleleAlt)
                        return count
                }
            }

            if(this.pos == (pos-alleleAlt.size()) && alleleType=="DEL" && obs=="-")
                return count
                
            ++count
        }
        return 0
    }
	
	Map toAnnovar(int alleleIndex=0) {
        Allele allele = this.getAlleles()[alleleIndex]
		String varType = allele.type
		if(varType == "DEL") {
			// For a deletion, Annovar returns "-"
            int altIndex = ref.lastIndexOf(allele.alt)+allele.alt.size()
			return [ pos: this.pos + altIndex, ref: this.ref.substring(altIndex), obs: "-" ]
		}
		else 
		if(varType == "INS") {
			return [ pos: this.pos + this.ref.size()-1, ref: "-", obs: allele.alt.substring(this.ref.size()) ]
		}
		else {
			return [ pos: this.pos, ref: this.ref, obs: this.alt ]	
		}
	}
    
    @CompileStatic
    Integer getTotalDepth(String sample) {
        int sampleIndex = this.header.samples.indexOf(sample)
        if(sampleIndex <0) // no genotype for sample?
            return 0
            
        if('DP' in this.header.formatMetaData) {
            return this.genoTypes[sampleIndex].DP as Integer
        }
        return (int) (0..alts.size()).collect { getAlleleDepths(it)[sampleIndex] }.sum()
    }
    
    /**
     * Return the number of reads supporting the given allele as a list
     * with one entry for each sample in the VCF
     * 
     * @param alleleIndex   index of allele, reference = 0
     */
    List<Integer> getAlleleDepths(int alleleIndex) {
        
        if('AD' in this.header.formatMetaData) {
            return this.genoTypes.collect { gt ->
                getGenotypeDepth(gt, alleleIndex)
            }
        }
        else
        // For a single sample VCF and a single alternate allele, we can use DP4 from samtools
        if(this.header.hasInfo('DP4') && alleleIndex<2 && this.header.samples.size() == 1) {
            List<Integer> dps = this.getInfo().DP4.tokenize(',')*.toInteger()
            
            assert dps.size() == 4 : "DP4 should always be list of 4 integers"
            
            int altDepth = dps[2] + dps[3]
            if(alleleIndex == 0) {
                return  [getInfo().DP.toInteger() - altDepth ] // arguably we should use the 1st and 2nd elements of DP4, but this
                                                          // seems to produce an unreliable result
            }
            else
            if(alleleIndex == 1) {
                return [altDepth]
            }
            else 
                assert false : "Only ref and single alternate allele supported!"
            
        }
        else {
            [0] * this.header.samples.size()
        }
    }
    
    @CompileStatic
    int getGenotypeDepth(Map<String,Object> gt, int alleleIndex) {
        List ad = (List)gt.AD
        if(!gt.AD) {
            if(alleleIndex==0 && (gt.GT=='0/0' || gt.GT=='0|0'))
                return (int)gt.DP
            else
                return 0i
        }
        else
        if(ad[alleleIndex]==null || (ad[alleleIndex] == "."))
            return 0i

        return (int)ad[alleleIndex]
    }
    
    /**
     * Return the distance of the balance of the default allele from 0.5
     * for the default alternate allele and first sample
     * 
     * @return
     */
    float getAlleleBalance() {
        return getAlleleBalance(1,0,0)
    }
    
    /**
     * Return the maximum allele balance for any sample for the default
     * alternate allele.
     * 
     * @return
     */
    float getMaxAlleleBalance() {
        (0..<header.samples.size()).collect { getAlleleBalance(1,0,it) }.max()
    }
    
    /**
     * Return the distance of the balance of the default allele from 0.5
     * for the specififed alternate allele (reference = 0)
     * 
     * @return
     */
    float getAlleleBalance(int allele1, int allele2, int sampleIndex=0) {
        
        // Allele balance only make sense for  heterozygous mutations: return 0
        if(this.getDosages(sampleIndex)[allele2] != 1) {
            return 0
        }
        
        List<Integer> ads1 = getAlleleDepths(allele1)
        List<Integer> ads2 = getAlleleDepths(allele2)
        
        if(ads1[sampleIndex] == 0)
            return 0
            
        if(ads2[sampleIndex] == 0)
            return 0
            
        return Math.abs(0.5 - ((float)ads1[sampleIndex] / (ads2[sampleIndex] + ads1[sampleIndex])))
    }
     
    /**
     * Return a JSON string representing the key details about this variant.
     * <b>Note:</b>
     * <li>If no sample is provided then the first sample is assumed
     * <li>The default information such as allele depths refers to the 
     *     first alternate allele. 
     * 
     * @return JSON string representing the variant
     */
    String toJson(String sample=null) {
        int sampleIndex = sample ? header.samples.indexOf(sample) : 0
        SnpEffInfo effect = getMaxEffect()
        Map effectInfo = [:]
        if(effect) {
            effectInfo = [ truncating : effect.isTruncating(), effect: effect.type, impact: effect.impact ]            
        }
        else {
            Map<String,Object> vepInfo = getMaxVep();
            if(vepInfo) {
                effectInfo = [
                    truncating : vepInfo.IMPACT == "HIGH",
                    effect : vepInfo.Consequence,
                    impact : vepInfo.IMPACT
                ]
            }
        }
        
        JsonOutput.toJson([
            chr : chr,
            alt : alt, 
            type: type,
            dosage: getDosages()[sampleIndex],
            alleles : this.getAlleles(),
            depths : [getAlleleDepths(0)[sampleIndex], getAlleleDepths(1)[sampleIndex]]
//            info : info
        ] + effectInfo)
    }
    
    /**
     * Return a list of each allele and its type.
     * Yes, VCF allows multiple types (INS,DEL) to be on the same line of 
     * a VCF file.
     * @deprecated  use getAlleles() instead
     */
    List<List<String>> getAllelesAndTypes() {
        return alts.collect {
            [it, convertType(ref,it)]
        }
    }
    
    /**
     * Return the index of the allele (if any) that matches the give other variant
     * @param other
     * @return
     */
    int findAlleleIndex(Allele other) {
        alleles.findIndexValues { me ->
            other.alt == me.alt && other.start == me.start && other.end == me.end && other.type == me.type
        }[0]
    }
    
    /**
     * Return a list of Allele objects representing alleles present on this line of the VCF
     * <p>
     * Note: currently this computes the results on the fly and they are not cached, so
     *       use with caution in computationaly intensive situations
     */
    List<Allele> getAlleles() {
        int index=0
        return alts.collect { obs ->
            String t = convertType(ref,obs)
            
            int start = pos
            int end = pos
            switch(t) {
                case "INS":
                    start = pos + ref.size()-1 // Insertion comes after reference sequence
                    end = start // since the positions are rel to reference genome, end must be same as start for an insertion
                    break
                case "DEL":
                    start = pos + obs.size()
                    end = pos + ref.size()
                    break
            }
            
            new Allele(index:index++, start: start, end: end, alt: obs, type: t)
        }
    } 
    
    /**
     * Return true if this variant overlaps the given range
     */
    boolean isCase(IRegion r) {
        return this.region.overlaps(r)
    }
    
    /**
     * Update the first alternate allele to the given value
     * 
     * @param alt   Alternate sequence of bases. Note that this should conform to VCF 
     *              spec, eg: for deletions it is expected at least one base of context 
     *              should be provided, and this should be consistent with the {@link #pos}
     *              field.
     */
    @CompileStatic
    void setAlt(String alt) {
        if(!alts) {
            alts = [alt]
        }
        else {
            alts[0] = alt
        }
        this.alt = alt
        type = convertType(ref,alt)
        
        if(!alt.isEmpty())
            altByte = (byte)alt.charAt(0)
        else
            altByte = '.'.charAt(0)
    }
    
    Region cachedRegion = null
    
    @CompileStatic
    Region getRegion() {
        if(!cachedRegion) {
            cachedRegion = new Region(this.chr, this.pos..(this.pos+this.size()))
        }
        return cachedRegion
    }
    
    /**
     * Return the list of genes impacted by this variant, in a manner that is
     * neutral to the annotator used (VEP and SnpEFF supported)
     * 
     * @param minVEPCons  the VEP impact level above which genes should
     *                      be returned
     */
    List<String> getGenes(String minVEPCons) {
        int thresholdConsequenceIndex = VEPConsequences.RANKED_CONSEQUENCES.indexOf(minVEPCons)
        if(this.header.getInfoMetaData("CSQ") || this.header.getInfoMetaData("ANN"))
            return getVepInfo().grep { vep ->
                AMPERSAND_SPLIT.split(vep.Consequence).every { VEPConsequences.RANKED_CONSEQUENCES.indexOf(it) <= thresholdConsequenceIndex }
            }*.SYMBOL
        else {
            return getSnpEffInfo().grep { snpeff ->
                String vepCons = SnpEffInfo.EFFECT_TO_VEP[snpeff.EFFECT_TO_VEP]
                if(!vepCons)
                    return false
                return (VEPConsequences.RANKED_CONSEQUENCES.indexOf(vepCons)) <= thresholdConsequenceIndex
            }*.gene
        }
    }
    
    @CompileStatic
    float getMaxVepMaf() {
       return (float)getVepInfo().collect { Map<String,Object> vep -> 
           ((List<String>)[vep.EA_MAF, vep.ASN_MAF, vep.EUR_MAF]).collect { String maf->
               maf?maf.split('&'):[]
           }.flatten().collect { Object maf -> ((String)maf).isFloat() ? ((String)maf).toFloat() : 0.0f }.max() ?: 0.0f }.max()
    }
    
    /**
     * Return the details of the most severe VEP consequence, as ranked by
     * {@link VEPConsequences#RANKED_CONSEQUENCES}.
     * 
     * @return  A map with key value pairs of VEP annotation fields.
     */
    Map<String,Object> getMaxVep() {
        def allVeps = getVepInfo().inject([]) { veps, vep ->
            veps.addAll(vep.Consequence.split("&").collect { [vep,it]})
            return veps;
        }
        return allVeps.max { 
            VEPConsequences.severityOf(it[1]) 
        }?.getAt(0)
    }
    
    /**
     * Return the consequence of the specified allele, at the moment,
     * from VEP annotations (later, from others). If multiple consequences
     * are present for the same allele, then the most severe consequence is 
     * returned.
     */
    String getConsequence(int alleleIndex) {
        def vep = getVepInfo()[alleleIndex]
        return vep?.Consequence?.split("&")?.max { VEPConsequences.severityOf(it)}
    }
    
    String getMaxVepImpact() {
        Map<String,Object> vep = getMaxVep()
        List maxCons = vep.Consequence?.tokenize('&')?.collect { VEPConsequences.VEP_IMPACTS[it]?:"UNKNOWN" }
        for(String impact in ["HIGH","MODERATE","MODIFIER","LOW"]) {
            if(impact in maxCons)
                return impact
        }
        return "UNKNOWN"
    }
    
}
