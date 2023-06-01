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

import java.text.NumberFormat
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

/**
 * This is deprecated and superseded by {@link VepConsequence}
 */
@Deprecated
class VEPConsequences {
    public static List<String> RANKED_CONSEQUENCES = [
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
 * Represents the diploid states of a variant.
 * Variants can only be in the `HEMI` state on the sex chromosomes.
 * 
 */
enum Genotype {
  Het,
  Hemi,
  Hom,
  Missing,
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
 * potentially multiple alleles and samples and corresponding annotations, genotypes, etc. 
 * However it is very common to work with normalised single-sample VCFs that guarantee only 
 * a single allele per site and represent just one sample.
 * To make operations on variants more convenient, some methods have both a sample/allele specific
 * and 'default allele' version. The default allele version will operate on the first allele
 * listed in the alleles for the first sample. These are convenient to use, but you should always keep in mind they
 * ignore any other alleles and samples and could be unsafe to use on non-normalised VCFs.
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
 * VCF vcf = VCF.parse("test.vcf")
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
    
    static final NumberFormat QUAL_FORMATTER = NumberFormat.getInstance()
    
    final static int FORMAT_FIELD_INDEX = 8
    final static int FIRST_GENOTYPE_INDEX = 9
    
    static {
        QUAL_FORMATTER.maximumFractionDigits = 2
        QUAL_FORMATTER.groupingUsed = false
    }
    
    /**
     * A specific allele in the context of a Variant in a VCF file
     * <p>
     * Each variant (line) in a VCF file may contain multiple alternate alleles.
     * This class represents one alternate allele in a VCF file.
     * 
     * @author simon.sadedin@mcri.edu.au
     */
    @CompileStatic
    final public class Allele {
        
        public Allele(int index, int start, int end, String alt, String type) {
            this.index = index;
            this.start = start;
            this.end = end;
            this.alt = alt;
            this.type = type;
        }

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
        
        String igv() {
            String url = "http://localhost:60151/load?locus=$chr:$start-$end"
            return new URL(url).text
        }
        
        String toString() { start != end ? "$start-$end $alt ($type)" : "$start $alt ($type)" }
    }
    
    
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
    
    /**
     * cached set of genotypes for the first sample.
     *
     */

    List<Genotype> sampleGenotypes

    Set<Pedigree> pedigrees
    
    List<Map<String,Object>> genoTypes
    
    /**
     * The original line from which this variant was parsed
     */
    String line
    
    List genoTypeFields
	
	
    @CompileStatic
	IntRange getRange() {
        if(isSV()) {
            return pos..pos+this.size()
        }
        else
    		return pos..<(pos+alts.max { it.size() }.size())
	}
    
    @CompileStatic
    private boolean parseFields() {
        parseFields(true)   
    }
    
    @CompileStatic
    Object asType(Class clazz) {
        if(clazz == Region)
            return new Region(this) 
        else
        if(clazz == Map) {
            Map infoMap = getInfo()
            return (Map)[
                chr: chr,
                pos: pos,
                ref: ref,
                alt: alt,
                qual: qual,
                dosages: this.getDosages()?.join(','),
                ad: this.getAlleleDepths(0)
            ] + getInfo() 
        }
    }
    
    /**
     * Parse the given line from a VCF file.
     * 
     * @param line
     */
    @CompileStatic
    private boolean parseFields(boolean ignoreHomRef) {
        
//        def fields = line.split('[\t ]{1,}')
        List<String>  fields = line.tokenize('\t').grep { it != null }
        if(ignoreHomRef && fields[4] == '<NON_REF>')
            return false
        
        chr = fields[0]
        id = fields[2]
        ref = fields[3]
        alts = COMMA_SPLIT.split(fields[4])
        alt = alts[0]
        qual = fields[5].equals(".") ? -1.0f : Float.parseFloat(fields[5])
        filter = fields[6]
        info = fields[7]
        
        // These lines are setting up default values to handle a "sites only" VCF without lots of special
        // downstream handling for the case where there is no genotype / sample
        boolean dirty = false
        if(fields[FORMAT_FIELD_INDEX] == null) { 
            fields[FORMAT_FIELD_INDEX] = 'GT'
            dirty = true
        }
            
        if(fields[9] == null)  {
            fields[9] = '0/1'
            dirty = true
        }
        
        if(dirty) {
            this.line = fields.join('\t')
        }
         
        pos = Integer.parseInt(fields[1])
        parseGenotypes(fields)
        
        setAlt(alt)
        return true
    }
    
    @CompileStatic
    private void parseGenotypes(List<String> fields) {
        if(fields.size() > FORMAT_FIELD_INDEX) {
            
          if(genoTypeFields == null)
              genoTypeFields = fields[FORMAT_FIELD_INDEX].tokenize(':')
          
          // Split the genoTypes field into separate values and parse them out
          if(fields.size()>FIRST_GENOTYPE_INDEX) {
              genoTypes = fields[9..-1].collect { String gt -> 
                  parseGenoTypeFields(gt)
              } 
          }
          else {
              genoTypes = []
          }
        }
        else {
            genoTypeFields = ['GT']
            genoTypes = (List<Map<String,Object>>) [ 
                [ GT: '0/1'] 
            ]
        }
    }
    
    final static Set numericGTFields = ["DP","GQ"] as Set
    final static Set numericListFields = ["AD"] as Set
    
    // Causes strange typecast error (Char => CharSequence at marked line below
    @CompileStatic
    Map<String,Object> parseGenoTypeFields(String gt) {
        [genoTypeFields, gt.tokenize(':')].transpose().collectEntries { Object fieldObj ->
              List field = (List)fieldObj
              String key = (String)field[0]
              String value = String.valueOf(field[1])
              if(key in numericGTFields) { 
                  return [key,convertNumericValue(value,null)] 
              }
              else
              if(key in numericListFields) {
                  return [key, value.equals('.') ? null : value.tokenize(",").collect { convertNumericValue(it,null) }]  // strange typecast exc under compilestatic
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
        
        if(value.isFloat())
            return value.toFloat()

        return defaultValue
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
    @CompileStatic
    void update(String desc, Closure c) {
        this.snpEffDirty = false
        if(line == null)
            line = [chr,pos,id?:".",ref,alts.join(","),".","PASS","ADDED"].join("\t")
            
        c(this)
        
        List<String> fields = TAB_SPLIT.split(line) as List
        fields[0] = chr
        fields[1] = String.valueOf(pos)
        fields[3] = ref
        fields[4] = alts?alts.join(','):"."
//        fields[5] = String.format("%2.2f",qual)
        fields[5] = QUAL_FORMATTER.format(qual)
        
        fields[6] = filter
        
        if(snpEffDirty) {
            // The only operation we support on the snpEff meta data is to remove some.
            // Hence we just need to rebuild it by concatenating it
            getInfo().EFF = this.snpEffInfo*.info.join(",")
        }
        
        fields[7] = getInfo()?.collect {k,v -> v != null?"$k=$v":k}.join(';')?:'.'
        
        if(this.@genoTypes != null) {
            rebuildGenotypes(fields)
        }
        
        fields[FORMAT_FIELD_INDEX] = genoTypeFields?.join(':')?:''
        
        line = fields.collect { it?:'' }.join('\t')
        
        // Check if we need to add an INFO header line
        if(this.header != null) {
            getInfo().collect { k,v ->
                if(!header.hasInfo(k)) {
                    header.addInfoHeader(k,desc,v)
                }
           }
        }
        this.cachedRegion = null
        this.cachedSize = null
        this.@info = null
    }
    
    /**
     * Return list of genotypes (state of variant) for 
     * each sample
     * @return
     */
    @CompileStatic
    List<Genotype> getSampleGenotypes(){
        if(sampleGenotypes != null)
            return sampleGenotypes

        List<String> gts = (List<String>)genoTypes*.GT
        List<Genotype> result = (List<Genotype>)gts.collect{
            if(it.isInteger()){
                return Genotype.Hemi
            }else{
                String[] split = PIPE_OR_SLASH_SPLIT.split(it)
                if(split[0] == split[1]){
                    if(split[0]=="."){
                        return Genotype.Missing
                    }
                    return Genotype.Hom
                }else{
                    return Genotype.Het
                }
            }
        }
        
        sampleGenotypes = result
            
        return result
    }

    /**
     * checks if the variant is heterozygous for the first sample
     * 
     * @return true if the variant call showed two distinct alleles
     */
    @CompileStatic
    boolean isHet() {
        if(this.getSampleGenotypes()[0] == Genotype.Het)
                return true
        return false
        }

    /**
     * tests the genotype of the variant in the specified sample
     *
     * @return true if the variant call showed two distinct alleles
     */
    @CompileStatic
    boolean isHet(int sampleIndex) {
        if(this.getSampleGenotypes()[sampleIndex] == Genotype.Het)
            return true
        return false
    }

            

    /**
     * tests the genotype of the variant in the specified sample
     *
     * @return true if the variant call showed two distinct alleles
     */
    @CompileStatic
    boolean isHet(String sampleName) {
        if(this.header == null)
            throw new IllegalStateException("Variant must have a header to query genotypes by sample name")
        int sampleIndex = header.samples.indexOf(sampleName)

        if(this.getSampleGenotypes()[sampleIndex] == Genotype.Het)
            return true
        return false
    }
    
    /**
     * tests the genotype of the variant in the first sample
     * 
     * @return true if variant call showed two identical alleles
     */
    @CompileStatic
    boolean isHom() {
        if(this.getSampleGenotypes()[0] == Genotype.Hom)
            return true
        return false
    }

    /**
     * tests the genotype of the variant in the specified sample
     * 
     * @return true if variant call showed two identical alleles
     */
    @CompileStatic
    boolean isHom(int sampleIndex) {
        if(this.getSampleGenotypes()[sampleIndex] == Genotype.Hom)
                return true
        return false
    }

    /**
     * tests the genotype of the variant in the specified sample
     *
     * @return true if the variant call showed two distinct alleles
     */
    @CompileStatic
    boolean isHom(String sampleName) {
        if(this.header == null)
            throw new IllegalStateException("Variant must have a header to query genotypes by sample name")
        int sampleIndex = header.samples.indexOf(sampleName)
        
        if(this.getSampleGenotypes()[sampleIndex] == Genotype.Hom)
            return true
        return false
    }


    /**
     * tests the genotype of the variant in the first sample
     *
     * @return true if variant call showed two identical alleles or only one allele
     */
    @CompileStatic
    boolean isHemiOrHom() {
        return this.isHet() || this.isHom()
    }

    /**
     * tests the genotype of the variant in the specified sample
     *
     * @return true if variant call showed two identical alleles or only one allele
     */
    @CompileStatic
    boolean isHemiOrHom(int sampleIndex) {
        return this.isHet(sampleIndex) || this.isHom(sampleIndex)
    } 

    /**
     * tests the genotype of the variant in the specified sample
     *
     * @return true if variant call showed two identical alleles or only one allele
     */
    @CompileStatic
    boolean isHemiOrHom(String sampleName) {
        if(this.header == null)
            throw new IllegalStateException("Variant must have a header to query genotypes by sample name")
        int sampleIndex = header.samples.indexOf(sampleName)

        return this.isHet(sampleIndex) || this.isHom(sampleIndex)
    } 

    /**
     * Update the per-sample genotype info fields by rebuilding them
     * from the genoTypes map.
     * 
     * @param fields    the parsed fields of the VCF line
     */
    @CompileStatic
    private void rebuildGenotypes(List fields) {
        assert header != null
        
        if(header.samples == null) {
            fields[9] = ''
            return
        }
        
        for(int i=0; i<this.header.samples.size(); ++i) {
            fields[9+i] = genoTypeFields.collect { fieldName ->
                def field = genoTypes[i][fieldName]
                if(field == null)
                     return '.'
                else
                if(field instanceof List)
                     return field.collect { it.is(null) ? '.' : it }.join(',') 
                 else
                     return field
            }.join(':')
        }
    }
    
    /**
     * VCF requires genotype to be the 1st field, but some tools (R, grrrr) write it in whatever order they 
     * feel like. This function corrects it.
     */
    void fixGenoTypeOrder() {
        def fields = line.split('[\t ]{1,}')
        
        def gtFields = fields[FORMAT_FIELD_INDEX]
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
            infos =  parseInfoString(this.info)        
        }
        return infos
    }
    
    /**
     * Parse the value from a VCF INFO field into a Map structure
     * <p>
     * Note this does not do a type safe parsing; all values are parsed to
     * Strings
     *  
     * @param value
     * @return
     */
    @CompileStatic
    static Map<String,Object> parseInfoString(final String value) {
        Map result = [:]
        if((value != null) && (value != '.')) {
            for(String s in value.tokenize(';')) {
                int i = s.indexOf('=')
                if(i<0) {
                    result[s]=Boolean.TRUE
                    continue
                }
                result[s.substring(0,i)] = s.substring(i+1)
            }
        }
        return result
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
   
    @CompileStatic
    List<Map<String,Object>> getVepInfo() {
        if(this.header == null)
            throw new IllegalStateException("Variant must have a header to query VEP information")
            
        String csqs = getInfo().CSQ 
        def vepFields = null
        if(csqs) {
           vepFields = this.header.getVepColumns("CSQ")
        }
        else
        if((csqs = getInfo().ANN)) {
           vepFields = this.header.getVepColumns("ANN")
        }
        else
        if((csqs = getInfo().vep)) {
           vepFields = this.header.getVepColumns("vep")
        }
        
        if(!csqs)
            return []
            
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
    @CompileStatic
    int getDosage() {
        getDosages(0)[0]
    }
    
    /**
     * Return list of dosages (number of copies of allele) for 
     * each sample for the first alternate allele.
     * @return
     */
    @CompileStatic
    List<Integer> getDosages() {
        getDosages(0)
    }
    
    /**
     * Return the number of copies of the given alternate allele for each sample in the VCF
     * 
     * <i>Note</i>: the first alternate allele is 0.
     */
    @CompileStatic
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
            
        List<String> gts = (List<String>)genoTypes*.GT
        List<Integer> result = (List<Integer>)gts.collect{PIPE_OR_SLASH_SPLIT.split(it)}*.count { 
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
            
        // If there are no samples then we assume it's a site-only VCF and return 1 for the dosage
        if(this.header.samples.isEmpty())
            return 1
            
        int sampleIndex = this.header.samples.indexOf(sampleName)
        if(sampleIndex < 0)
            throw new IllegalArgumentException("Sample $sampleName not found in VCF. Known samples are $header.samples")
            
        List<Integer> allDosages = getDosages()
        if(allDosages.isEmpty()) // Another workaround for sites-only VCF
            return 1
        
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
            
        String prefix = sequence.substring(0,position-1)
        String suffix = sequence.substring(position -1 + ref.size())
        
//        println "Ref   : " + ref
//        println "Sequen: " + sequence
//        println "Prefix: $prefix"
//        println "Introd: $alt"
//        println "Suffix: $suffix"
//        
        prefix + alt + suffix
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
//        return (alt == "DEL" || alt == "<DEL>" || alt == "DUP" || alt == "<DUP>" || alt == "INV" || alt == "<INV>") 
        return (alt == "DEL" || alt == "DUP" || alt == "INV" || 
            (alt.startsWith('<') && alt.endsWith('>')) || 
            (alt.size()>1 && (alt[1] == '[' || alt[-2] == ']'))
        )
    }
    
    Integer cachedSize = null
    
    /**
     * Returns the change in size of the genome caused by this variant's default allele.
     * <p>
     * Note: since point mutations / SNVs don't change the size, they have size zero
     * 
     * @return difference in genome size caused by this variant's default allele
     */
    @CompileStatic
    int size() {
        
        if(!cachedSize.is(null))
            return cachedSize
        
        if(isSV()) {
           Object svLen = this.getInfo().SVLEN
           if(svLen) {
                cachedSize = Math.abs(Integer.parseInt(String.valueOf(svLen)))
                return cachedSize
           }
                   
           String endValue = (String)this.getInfo().END
           String chr2 = this.getInfo().CHR2

           // Can only use END if on the same chromosome
           if((chr2.is(null) || (chr2 == this.chr)) && endValue) {
               int end = Integer.parseInt(endValue)
               cachedSize = end - pos
               return cachedSize
           }
           return 0
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
        if(!parsed.parseFields())
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
    Integer getTotalDepth() {
        if('DP' in this.header.formatMetaData) {
            return (this.genoTypes[0].DP?:0) as Integer
        }
        return (int) (0..alts.size()).collect { getAlleleDepths(it)[0] }.sum()
    } 
    
    @CompileStatic
    Integer getTotalDepth(String sample) {
        int sampleIndex = this.header.samples.indexOf(sample)
        if(sampleIndex <0) // no genotype for sample?
            return 0
            
        return this.getTotalDepth(sampleIndex)
    }
    
    @CompileStatic
    Integer getTotalDepth(final int sampleIndex) {
        if('DP' in this.header.formatMetaData) {
            return (this.genoTypes[sampleIndex].DP?:0) as Integer
        }
        return (int) (0..alts.size()).collect { getAlleleDepths(it)[sampleIndex] }.sum()
    }
    
    /**
     * Return the depth of the first alternate allele for the first sample
     * @return
     */
    @CompileStatic
    int getAltDepth() {
        return getAlleleDepths(1)[0]
    }
    
    /**
     * Return the number of reads supporting the given allele as a list
     * with one entry for each sample in the VCF
     * 
     * @param alleleIndex   index of allele, reference = 0
     */
    @CompileStatic
    List<Integer> getAlleleDepths(int alleleIndex) {
        
        boolean hasAD = ('AD' in this.header?.formatMetaData) ||
                        ('AD' in this.genoTypes[0])
        
        if(hasAD) {
            return this.genoTypes.collect { gt ->
                getGenotypeDepth(gt, alleleIndex)
            }
        }
        else
        // For a single sample VCF and a single alternate allele, we can use DP4 from samtools
        if(this.header?.hasInfo('DP4') && alleleIndex<2 && this.header?.samples.size() == 1) {
            List<Integer> dps = ((String)this.getInfo()['DP4']).tokenize(',')*.toInteger()
            
            assert dps.size() == 4 : "DP4 should always be list of 4 integers"
            
            int altDepth = dps[2] + dps[3]
            if(alleleIndex == 0) {
                return  [((String)getInfo()['DP']).toInteger() - altDepth ] // arguably we should use the 1st and 2nd elements of DP4, but this
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
    
    /**
     * Tries to intelligently parse the genotypes in this variant by converting all
     * fields that look numeric to a numeric form.
     * <p>
     * NOTE: this is not done respecting information in the VCF header, and may return
     * inconsistent data types for different variants for the same field - eg: if a number
     * appears to be an integer, it will ve converted to an integer, even if the 
     * field may be sometimes fractional.
     */
    List<Map> getParsedGenotypes() {
        this.genoTypes.collect { Map gtInfo ->
            gtInfo.collectEntries { k,x ->
                def value = x
                if(x instanceof String)  {
                    def listValue = Utils.toNumberList(x)
                    if(listValue)
                        value = listValue
                    else
                        value = Utils.convertNum(x)
                }

                [k, value]
            }
        }
    }
    
    Map<String,Object> getParsedInfo() {
        this.getInfo().collectEntries { k,x ->
            def value = x
            if(x instanceof String)  {
                def listValue = Utils.toNumberList(x)
                if(listValue)
                    value = listValue
                else
                    value = Utils.convertNum(x)
            }

            [k, value]
        }
    }
    
    @CompileStatic
    int getGenotypeDepth(Map<String,Object> gt, int alleleIndex) {
        List ad = (List)gt.AD
        if(!gt.AD) {
            if(alleleIndex==0 && (gt.GT=='0/0' || gt.GT=='0|0'))
                return (int)(gt.DP?:0i)
            else
                return 0i
        }
        else
        if(ad[alleleIndex]==null || (ad[alleleIndex].equals(".")))
            return 0i

        return (int)ad[alleleIndex]
    }
    
    /**
     * Return the distance of the balance of the default allele from 0.5
     * for the default alternate allele and first sample
     * 
     * @return
     */
    @CompileStatic
    float getAlleleBalance() {
        return getAlleleBalance(1,0,0)
    }
    
    /**
     * Return the maximum allele balance for any sample for the default
     * alternate allele.
     * 
     * @return
     */
    @CompileStatic
    float getMaxAlleleBalance() {
        (0..<header.samples.size()).collect { getAlleleBalance(1,0,it) }.max()
    }
    
    /**
     * Return the distance of the balance of the default allele from 0.5
     * for the specififed alternate allele (reference = 0)
     * 
     * @return
     */
    @CompileStatic
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
            
        return Math.abs(0.5f - ((float)ads1[sampleIndex] / (ads2[sampleIndex] + ads1[sampleIndex])))
    }
    
    /**
     * Return the fraction of reads supporting the first alternate alelle for the
     * first sample in the VCF
     */
    @CompileStatic
    float getVaf() {
        return this.getVaf(1)
    }
   
    /**
     * Return the fraction of reads supporting the given allele for the specified 
     * sample
     * <p>
     * The 0th allele is the reference, so typically you would want to use
     * 1 or more.
     * 
     * @param alleleIndex   index of the allele to return the frac
     * @param sampleIndex   index of sample (in order of VCF header)
     * @return  value between zero and 1 indicating fraction of reads supporting allele
     */
    @CompileStatic
    float getVaf(int alleleIndex, int sampleIndex=0) {
        List<Integer> ads = getAlleleDepths(alleleIndex)
        final int totalDepth = this.getTotalDepth(sampleIndex)
        if(totalDepth == 0)
            return 0f
        return ads[sampleIndex] / totalDepth
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
    @CompileStatic
    int findAlleleIndex(Allele other) {
        return this.alleles.findIndexOf { Allele me ->
            other.alt == me.alt && other.start == me.start && other.end == me.end && other.type == me.type
        }
    }
    
    /**
     * Return a list of Allele objects representing alleles present on this line of the VCF
     * <p>
     * Note: currently this computes the results on the fly and they are not cached, so
     *       use with caution in computationaly intensive situations
     */
    @CompileStatic
    List<Allele> getAlleles() {
        int index=0
        return alts.collect { String obs ->
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
            
            new Allele(index, start, end, obs, t)
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
        int thresholdConsequenceIndex = VepConsequence.fromTerm(minVEPCons).ordinal()
        if(this.header.getInfoMetaData("CSQ") || this.header.getInfoMetaData("ANN"))
            return getVepInfo().grep { vep ->
                AMPERSAND_SPLIT.split(vep.Consequence).every { VepConsequence.fromTerm(it).ordinal() <= thresholdConsequenceIndex }
            }*.SYMBOL
        else {
            return getSnpEffInfo().grep { snpeff ->
                String vepCons = SnpEffInfo.EFFECT_TO_VEP[snpeff.EFFECT_TO_VEP]
                if(!vepCons)
                    return false
                return (VepConsequence.fromTerm(vepCons).ordinal()) <= thresholdConsequenceIndex
            }*.gene
        }
    }
    
    @CompileStatic
    float getMaxVepMaf() {
       def result = getVepInfo().collect { Map<String,Object> vep -> 
           
           if(vep.containsKey('MAX_AF')) {
               return ((String)vep.MAX_AF).isFloat() ? ((String)vep.MAX_AF).toFloat() : 0.0f
           }
           
           ((List<String>)[vep.EA_AF, vep.EAS_AF, vep.AFR_AF, vep.AMR_AF, vep.EUR_AF]).collect { String maf->
               maf?maf.split('&'):[]
       }.flatten().collect { Object maf -> ((String)maf).isFloat() ? ((String)maf).toFloat() : 0.0f }.max() ?: 0.0f }.max()
           
       return result ? (float)result : 0.0f
    }
    
    /**
     * Return the details of the most severe VEP consequence, as ranked by
     * {@link VepConsequence} ordinal value.
     * 
     * @return  A map with key value pairs of VEP annotation fields.
     */
    @CompileStatic
    Map<String,Object> getMaxVep() {
        List<List> allVeps = (List<List>)getVepInfo().inject([]) { veps, vep ->
            String cons = vep.Consequence
            veps.addAll(cons.split("&").collect { [vep,it]})
            return veps;
        }
        return allVeps.max { List vepAndCons ->
            VepConsequence.severityOf(vepAndCons[1] as String)
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
        return vep?.Consequence?.split("&")?.max { VepConsequence.severityOf(it)}
    }
    
    @CompileStatic
    long getXpos() {
        XPos.computePos(this.chr, this.pos)
    }
    
    String getMaxVepImpact() {
        Map<String,Object> vep = getMaxVep()
        List maxCons = vep.Consequence?.tokenize('&')?.collect { VepConsequence.fromTerm(it).impact.name() }
        for(String impact in ["HIGH","MODERATE","MODIFIER","LOW"]) {
            if(impact in maxCons)
                return impact
        }
        return "UNKNOWN"
    }
    
    String igv() {
        alleles[0].igv()
    }
    
}
