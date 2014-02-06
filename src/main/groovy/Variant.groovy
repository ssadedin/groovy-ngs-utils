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

import groovy.transform.CompileStatic

/**
 * Support for parsing annotations from SnpEff
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class SnpEffInfo {
    String type
    String impact
    String gene
    String transcript
    
    // The original subclause within the info line from which
    // this SnpEffInfo object was parsed
    String info
    
    String toString() {
        "type=$type,gene=$gene,rank=$impact"
    }
}

/**
 * Help with parsing VCFs
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class Variant {
    
    String chr
    int pos
    String ref
    String alt
    String [] alts
    String id
    String type
    String info
    String filter
    float qual
    byte altByte
    VCF header

    Map<String,String> infos
    
    List<SnpEffInfo> snpEffInfo
    
    List<Integer> dosages
    
    Set<Pedigree> pedigrees
    
    List<Map<String,Object>> genoTypes
    
    String line
    
    String [] genoTypeFields
    
    // @CompileStatic
    void parseFields(String line) {
        
        def fields = line.split('[\t ]{1,}')
        
        chr = fields[0]
        id = fields[2]
        ref = fields[3]
        alts = fields[4].split(",")
        alt = alts[0]
        qual = (fields[5] == ".") ? -1.0f : Float.parseFloat(fields[5])
        filter = fields[6]
        info = fields[7]
        
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
                      return [it[0],it[1].split(",")*.toFloat()] 
                  }
                  else
                      it
              }
          } 
        }
          
        pos = Integer.parseInt(fields[1])
        
        type = convertType(ref,alt)
        altByte = (byte)alt.charAt(0)
    }
    
    /**
     * Convert reference and alternate allele strings into a 
     * mutation type, being one of "SNP","INS","DEL"
     * 
     * @param refSeq
     * @param altSeq
     * @return
     */
    String convertType(String refSeq, String altSeq) {
        String result = "SNP"
        if(refSeq.size() < altSeq.size())
            result = 'INS'
        else
        if(refSeq.size() > altSeq.size())
            result = 'DEL'
            
        return result
    }
    
    boolean snpEffDirty = false
    
    void update(Closure c) {
        update("Groovy Variant Processing " + (new Date()).toString(),c)
    }
    
    /**
     * Allows various fields to be updated and then synchronises the 
     * rest of the data with those updated fields
     * Not all fields are supported! see the ones that are set below.
     * <p>
     * The only update to snpEFF information is to remove individual 
     * annotations.
     */
    void update(String desc, Closure c) {
        this.snpEffDirty = false
        c()
        def fields = line.split('[\t ]{1,}')
        fields[0] = chr
        fields[1] = String.valueOf(pos)
        fields[3] = ref
        fields[4] = alt
        fields[5] = String.format("%2.2f",qual)
        
        if(snpEffDirty) {
            // The only operation we support on the snpEff meta data is to remove some.
            // Hence we just need to rebuild it by concatenating it
            getInfo().EFF = this.snpEffInfo*.info.join(",")
        }
        
        fields[7] = getInfo().collect { k,v -> "$k=$v"}.join(';')
        
        line = fields.join('\t')
        
        // Check if we need to add an INFO header line
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
            for(String s in info.split(';')) {
                int i = s.indexOf('=')
                if(i<0) {
                    infos[s]=Boolean.TRUE
                    continue
                }
                infos[s.substring(0,i)] = s.substring(i+1)
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
        def vepFields = this.header.vepColumns
        this.getInfo().CSQ.split(",").collect { csq -> [vepFields,csq.split("\\|")].transpose().collectEntries() }
    }
    
    /**
     * Return the most impactful SnpEff effect for this variant
     */
    SnpEffInfo getMaxEffect() {
        def info = getSnpEffInfo()
        // TODO: at the moment the max effect is just abitrarily chosen, but we know we could / should 
        // do better
        return info?.max { it.rank }
    }
    
    List<Integer> getDosages() {
        if(dosages != null)
            return dosages
        
        // TODO: dosages work differently in de novo callers, 
        // they will just supply the base instead of allele 
        // Here is how it was computed in another script    
        /*
          def alleles = s.split(':')[0].split('/')
          if(alleles[0] in ['A','C','T','G']) 
              alleles[0] == alleles[1] ? 2 : 1
          else
              alleles.collect { 
                  it == '.' ? -1 : Integer.parseInt(it) 
              }.countBy { it }[1]?:0
 
         */
        dosages = genoTypes*.GT*.split('/')*.sum { it == '.' ? 0 : Integer.parseInt(it) } 
    }

    @CompileStatic 
    int sampleDosage(String sampleName) {
        return getDosages()[header.samples.indexOf(sampleName)]
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
        return genoTypes[header.samples.indexOf(sampleName)]
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
        if(this.@pedigrees == null) {
          def dsgs = getDosages()
          this.@pedigrees = new HashSet()
          for(int i=0; i<dsgs.size(); ++i) {
              if(dsgs[i] > 0) {
                  def ped = header.findPedigreeBySampleIndex(i)   
                  if(!ped) 
                      System.err.println("WARNING: no pedigree information found for sample " + header.samples[i])
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
    
    int size() {
        Math.abs(ref.size() - alt.size())
    }

    static Variant parse(String line) {
        Variant parsed = new Variant(line:line)
        parsed.parseFields(line)
        return parsed
    }
    
    String toString() {
        "$chr:$pos $ref/$alt"
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
            
            def alleleAlt = alleleTypePair[0]
            def alleleType = alleleTypePair[1]
            
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
    
    /**
     * Return a list of each allele and its type.
     * Yes, VCF allows multiple types (INS,DEL) to be on the same line of 
     * a VCF file!
     */
    List<List<String>> getAllelesAndTypes() {
        return alts.collect {
            [it, convertType(ref,it)]
        }
    }
}
