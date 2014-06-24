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

import com.sun.xml.internal.ws.developer.MemberSubmissionAddressing;

import groovy.json.JsonOutput;
import groovy.transform.CompileStatic

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
        'NON_SYNONYMOUS_CODING',
        'EXON',
        'START_GAINED',
        'SYNONYMOUS_CODING',
        'SPLICE_SITE_ACCEPTOR',
        'SPLICE_SITE_DONOR',
        'INTRON',
        '3_PRIME',
        '5_PRIME',
        'DOWNSTREAM',
        'UPSTREAM',
        'INTRAGENIC'
    ]

    /**
     * One of INS, DEL or SNP
     */
    String type
    
    /**
     * The kind of impact. One of the values from {@link SnpEffInfo#EFFECTS_RANKED}
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
    
    String toString() {
        "type=$type,gene=$gene,rank=$impact"
    }
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
     * INS, DEL or SNP
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
 * @author simon.sadedin@mcri.edu.au
 */
class Variant implements IRegion {
    
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
     * The VCF to which this Variant is linked. This enables the variant to intelligently
     * parse INFO fields and be aware of sample names. If the header is not provided,
     * many functions still work, but some functions will be disabled.
     */
    VCF header

    Map<String,String> infos
    
    List<SnpEffInfo> snpEffInfo
    
    List<Integer> dosages
    
    Set<Pedigree> pedigrees
    
    List<Map<String,Object>> genoTypes
    
    /**
     * The original line from which this variant was parsed
     */
    String line
    
    String [] genoTypeFields
	
	
	IntRange getRange() {
        if(isSV()) {
            return pos..pos+this.size()
        }
        else
    		return pos..(pos+alts.max { it.size() }.size())
	}
    
    // @CompileStatic
    /**
     * Parse the given line from a VCF file.
     * 
     * @param line
     */
    private void parseFields(String line) {
        
        def fields = line.split('[\t ]{1,}')
        
        chr = fields[0]
        id = fields[2]
        ref = fields[3]
        alts = fields[4].split(",")
        alt = alts[0]
        qual = (fields[5] == ".") ? -1.0f : Float.parseFloat(fields[5])
        filter = fields[6]
        info = fields[7]
        
        pos = Integer.parseInt(fields[1])
        
        parseGenotypes(fields)
        
        setAlt(alt)
    }
    
    private void parseGenotypes(String [] fields) {
        if(fields.length > 8) {
            
          if(genoTypeFields == null)
              genoTypeFields = fields[8].split(':')
          
          // Split the genoTypes field into separate values and parse them out
          def numericGTFields = ["DP","GQ"] as Set
          def numericListFields = ["AD"] as Set
          genoTypes = fields[9..-1].collect { String gt -> 
              [genoTypeFields, gt.split(':')].transpose().collectEntries { 
                  if(it[0] in numericGTFields) 
                      [it[0],it[1].toFloat()] 
                  else
                  if(it[0] in numericListFields) {
                      return [it[0],it[1].split(",")*.replaceAll("\\.","0")*.toFloat()] 
                  }
                  else
                      it
              }
          } 
        }
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
            
        c()
        
        def fields = line.split('\t')
        fields[0] = chr
        fields[1] = String.valueOf(pos)
        fields[3] = ref
        fields[4] = alt?:"."
        fields[5] = String.format("%2.2f",qual)
        
        if(snpEffDirty) {
            // The only operation we support on the snpEff meta data is to remove some.
            // Hence we just need to rebuild it by concatenating it
            getInfo().EFF = this.snpEffInfo*.info.join(",")
        }
        
        fields[7] = getInfo().collect {k,v -> v != null?"$k=$v":k}.join(';')
        
        line = fields.join('\t')
        
        // Check if we need to add an INFO header line
        if(this.header != null) {
            getInfo().collect { k,v ->
                
                if(!header.headerLines.find { it.startsWith("##INFO=<ID=$k") }) {
                    
                    int lastInfo = header.headerLines.findLastIndexOf { it.startsWith("##INFO=") }
                    if(lastInfo < 0)
                        lastInfo = 1
                    
                    String valueType = "String"
                    if(v instanceof Integer) {
                        valueType = "Integer"
                    }
                    else
                    if(v instanceof Float) {
                        valueType = "Float"
                    }
                    else
                    if(v instanceof Double) {
                        valueType = "Double"
                    }
                     
                    header.headerLines = header.headerLines[0..lastInfo] + 
                        ["##INFO=<ID=$k,Number=1,Type=${valueType},Description=\"$desc\">"] +
                        header.headerLines[(lastInfo+1)..-1]
                }
           }
        }
        this.@info = null
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
    Map<String,String> getInfo() {
        if(infos == null) {
            infos =  [:]
            if(info != null) {
                for(String s in info.split(';')) {
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
            
        def genes = (eff =~ /[A-Z]*\(([^\)]*)\)/).collect { (it[1]+" ").split(/\|/)[5] }
        def txs = (eff =~ /[A-Z]*\(([^\)]*)\)/).collect { (it[1]+" ").split(/\|/)[8] }
        def effs = (eff =~ /([A-Z0-9_]*)\(([^\)]*)\)/).collect { it[1] }
        def ranks = (eff =~ /[A-Z]*\(([^\)]*)\)/).collect { (it[1]+" ").split(/\|/)[0] }
        def infos = eff.split(',')
        
        snpEffInfo = [genes,effs,ranks,infos,txs].transpose().collect { new SnpEffInfo(gene:it[0], type: it[1], impact:it[2], info: it[3], transcript:it[4]) }
    }
    
    List<Map<String,Object>> getVepInfo() {
        if(this.header == null)
            throw new IllegalStateException("Variant must have a header to query VEP information")
            
        def vepFields = this.header.vepColumns
        this.getInfo().CSQ.split(",").collect { csq -> [vepFields,csq.split("\\|")].transpose().collectEntries() }
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
        def result = genoTypes*.GT*.split('/')*.count { 
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
        int sampleIndex = header.samples.indexOf(sampleName)
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
        return (alt == "DEL" || alt == "<DEL>" || alt == "DUP" || alt == "<DUP>") 
    }
    
    @CompileStatic
    int size() {
        if(alt == "DEL" || alt == "<DEL>") {
           Object svLen = this.getInfo().SVLEN
           if(!svLen)
               throw new RuntimeException("VCF file contains structural variants but does not have SVLEN information in INFO field")
           
            return Math.abs(svLen.toInteger())
        }
        else
        if(alt == "DUP" || alt == "<DUP>") {
           Object svLen = this.getInfo().SVLEN
           if(!svLen)
               throw new RuntimeException("VCF file contains structural variants but does not have SVLEN information in INFO field")
            return svLen.toInteger()            
        }
        else
        return Math.abs(ref.size() - alt.size())
    }

    static Variant parse(String line) {
        Variant parsed = new Variant(line:line)
        parsed.parseFields(line)
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
              
            if(this.pos == (pos-this.ref.size()+1) && alleleType=="INS" && alleleAlt.endsWith(obs))
                return count
                
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
    
    String toJson() {
        JsonOutput.toJson([
            chr : chr,
            alt : alt, 
            type: type,
            dosage: getDosages()[0],
            alleles : this.getAlleles(),
            info : info
        ])
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
            
            new Allele(start: start, end: end, alt: obs, type: t)
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
    void setAlt(String alt) {
        if(!alts) {
            alts = [alt]
        }
        else {
            alts[0] = alt
        }
        this.alt = alt
        type = convertType(ref,alt)
        altByte = (byte)alt.charAt(0)
    }
    
    Region cachedRegion = null
    
    @CompileStatic
    Region getRegion() {
        if(!cachedRegion) {
            cachedRegion = new Region(this.chr, this.pos..(this.pos+this.size()))
        }
        return cachedRegion
    }
}
