/*
 *  Groovy NGS Utils - Some simple utilites for processing Next Generation Sequencing data.
 *
 *  Copyright (C) 2014 Simon Sadedin, ssadedin<at>gmail.com
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
package gngs.db 

import groovy.sql.Sql
import groovy.util.logging.Log;

import com.xlson.groovycsv.PropertyMapper;

import gngs.db.Schema;
import gngs.*
import gngs.Variant.Allele

/**
 * VariantDB implements a simple, embedded, variant tracking database that
 * can be used to track variants identified in multiple samples and 
 * sequencing runs. It allows quick computation of variant counts
 * across samples, families, batches of samples, etc.
 * <p>
 * For the full power of the database, you shoudl use a PED file, parsed
 * with the Pedigrees class when adding variants and samples. However this
 * parameter can be passed as null, in which case all samples are treated
 * as singletons. You should be aware that variant counts will not be family
 * aware in such cases and thus will be distorted if your sequencing has
 * large pedigrees and compared to the overall sample count.
 * <p>
 * To add variants, parse them using the VCF class and then simple use the 
 * #add method:
 * <pre>
 * VariantDB db = new VariantDB("test.db")
 * VCF.parse("test.vcf") { v ->
 *     db.add("batch1", null, v)
 * }
 * </pre>
 * 
 * 
 * @author simon
 */
@Log
class VariantDB {
    
    /**
     * The file name of the database that is connected to
     */
    String connectString
    
    /**
     * The file name of the database that is connected to
     */
    String driver
    
    /**
     * The actual database connection
     */
    Sql db
    
    /**
     * Schema for database
     */
    Schema schema 
    
    final static ANNOVAR_ANNOTATION_PROFILE = [
        ONEKG_FREQ : ["1000g2010nov_all", "1000g2010nov_ALL"],
        ESP_FREQ : ["ESP5400_all", "ESP5400_ALL"],
        EXAC_FREQ : ["exac03"],
        AACHANGE : ["AAChange", "AAChange_RefSeq", "AAChange.RefSeq"],
        DBSNP : ["dbSNP138", "dbsnp138"]
    ]
    
    Map<String, List<String>> annotationProfile = ANNOVAR_ANNOTATION_PROFILE
    
    /**
     * Open or create a VariantDB using default settings and the given
     * file name. These connect to a SQLite database with the given name.
     * 
     * @param fileName
     */
    public VariantDB(String fileName) {
        this.connectString = "jdbc:sqlite:$fileName"
        this.driver = "org.sqlite.JDBC"
        this.schema = new Schema()
        this.init()
    }
    
    /**
     * Check the database exists and upgrade it if necessary
     */
    void init() {
        
        log.info "Creating database connection to $connectString"
        Class.forName(this.driver)
        this.db = Sql.newInstance(connectString)
        
        // For performance
        db.execute("PRAGMA synchronous = 0;")
        db.execute("PRAGMA journal_mode = WAL;")
        this.schema.checkSchema(db)
    }

    /**
     * Return row of database representing sample, or null if it does not
     * exist.
     */
    def findSample(String sampleId) {
        sampleId = trimSampleId(sampleId)
        db.firstRow("select * from sample where sample_id = $sampleId")
    }
    
    String trimSampleId(String sampleId) {
        sampleId.replaceAll('_S[0-9]*$','')
    }
    
    /**
     * Add information about the given sample to the database
     * 
     * @param sampleId
     * @param peds
     * @return
     */
    def addSample(String sampleId, Pedigrees peds=null, String cohort=null, String batch=null) {
        Pedigree family = peds ? peds.subjects[sampleId] : null
        Subject subject = family?.individuals?.find { it.id == sampleId }
        if(family && !subject) 
            throw new IllegalStateException("Sample $sampleId could not be located in the pedigree file. Known samples are ${peds.subjects*.key}")
            
        sampleId = trimSampleId(sampleId)
                
        db.execute("""
            insert into sample (id, sample_id, father_id, mother_id, family_id, phenotype, cohort, batch, created) 
                        values (NULL, 
                                $sampleId, 
                                ${family?.motherOf(sampleId)?.id}, 
                                ${family?.fatherOf(sampleId)?.id}, 
                                ${family?.id?:sampleId}, 
                                ${subject?.phenoTypes?.getAt(0)?:0}, 
                                ${cohort},
                                ${batch},
                                datetime('now'));
        """)
        return findSample(sampleId)
    }
    
    /**
     * Find a variant in the database by its start position
     */
    def findVariant(Variant variant, Allele allele) {
        db.firstRow("select * from variant where chr=$variant.chr and start=$allele.start and alt=$allele.alt")
    }
    
    /**
     * Returns a map containing the following keys:
     *   <li>sampleCount - the number of unique samples in which the variant was observed
     *   <li>familyCount - the number of unique families in which the variant was observed
     *   
     * @param chr   Chromosome of variant
     * @param start starting position of DNA change caused by variant (note this may be different 
     *        to VCF position depending on the representation of your indel in the VCF!)
     * @param end   end position of DNA change caused by variant
     * @param alt   alternate sequence
     * @return  Map with keys indicating counts of observations of this variant
     */
    Map countObservations(String chr, int start, int end, String alt) {
        Variant v = new Variant(chr: chr, ref:"N", alt:alt, pos: start)
        countObservations(v,v.alleles[0])
    }
    
    /**
     * Return a counts of the number of observations of the given variant.
     * Two values are returned: samples and families. 
     */
    Map countObservations(Variant v, Allele allele) {
        int sampleCount = db.firstRow(
           """select count(distinct(o.sample_id))
                       from variant_observation o,
                            variant v 
                       where o.variant_id = v.id
                         and v.chr = $v.chr
                         and v.start = $allele.start
                         and v.end = $allele.end
                         and v.alt = $allele.alt
        """)[0]
        
        def familySql = """select count(distinct(s.family_id))
                     from variant_observation o,
                          variant v, 
                          sample s
                     where o.variant_id = v.id
                         and v.chr = $v.chr
                         and v.start = $allele.start
                         and v.end = $allele.end
                         and v.alt = $allele.alt
                         and o.sample_id = s.id
        """
        
        println(familySql)
        int familyCount = db.firstRow(familySql)[0]
        
        return [samples: sampleCount, families: familyCount]
    }
    
    /**
     * Search for the given value in the annotations as a frequency
     * value.
     * 
     * @param type          The type of value to search for, must be one of the 
     *                      predefined annotation profile keys
     * @param annotations   The annotations to search
     * @return  the value of the annotation or null if it is not available
     */
     Object getAnnotation(String type, annotations) {
       def key = annotationProfile[type].find { key -> try {annotations[key]} catch(Exception e) { null } }
       if(key)
           return annotations[key]
       else
           return null
     }
    
    /**
     * Search for the given value in the annotations as a frequency
     * value.
     * 
     * @param type          The type of frequency to search for, must be a key in the annotation profile
     * @param annotations   The annotations to search
     * @return
     */
    Float getFreq(String type, annotations) {
        def value = getAnnotation(type,annotations)
        if(!value)
            return 0.0f
        if(value instanceof Float)
            return value
        if(value == ".") 
            return 0.0f
        return value.isFloat() ? value.toFloat() : 0.0f
    }
    
    /**
     * Add the given variant to the database. If sampleToAdd is specified,
     * add it for only the given sample. Otherwise add it for all the samples that 
     * are genotyped to have the variant
     * 
     * @param annotations   optional annotations to draw from. These can be used when 
     *                      the annotations are not embedded in the VCF file. The annotations
     *                      are queried for keys defined in the annotationProfile field.
     *                      The default mappings are set up to look for Annovar annotations.
     *                      The annotations themselves can be any object having properties which 
     *                      will be queried. In practise, a CSV parser PropertyMapper is
     *                      what is being used here to pass in Annovar annotations.
     * 
     * @return count of variants added
     */
    int add(String batch, Pedigrees peds, Variant v, Allele alleleToAdd=null, String cohort=null, String sampleToAdd=null, def annotations=null) {
        
        int countAdded = 0
        def alleles = alleleToAdd ? [alleleToAdd] : v.alleles
        for(Allele allele in v.alleles) {
            def variant_row = findVariant(v, allele)
            if(!variant_row) {
                
                if(annotations == null) {
                    def vep = v.getVepInfo()[allele.index]
                    def cons = v.getConsequence(allele.index)
                    if(cons != null && cons.indexOf("(")>=0) 
                        cons = cons.substring(0, cons.indexOf("("))
                    db.execute("""
                        insert into variant (id,chr, pos, start,end,ref,alt, sift, polyphen, condel, consequence, max_freq,dbsnp_id) 
                                    values (NULL, $v.chr, $v.pos, $allele.start, $allele.end, ${v.ref}, $allele.alt, ${vep?.SIFT}, ${vep?.PolyPhen}, ${vep?.Condel}, 
                                           ${cons}, $v.maxVepMaf, $v.id);
                    """)
                }
                else {
                    
                    float maxFreq = ["ONEKG_FREQ", "ESP_FREQ", "EXAC_FREQ"].collect { getFreq(it,annotations) }.max()
                    String aaChange = getAnnotation("AACHANGE", annotations)
                    db.execute("""insert into variant (id,chr,pos,start,end,ref,alt,consequence,protein_change,max_freq, dbsnp_id) 
                                   values (NULL, $v.chr, 
                                                 $v.pos, 
                                                 $allele.start, $allele.end, 
                                                 $v.ref, $allele.alt,
                                                 $annotations.ExonicFunc,
                                                 $aaChange,
                                                 $maxFreq, 
                                                 ${getAnnotation("DBSNP", annotations)})
                               """)
                }
            }
            variant_row = findVariant(v,allele)
            
            // For every sample carrying the allele, add an observation
            for(String sampleId in v.header.samples.grep { v.sampleDosage(it) }) {
                
                if(sampleToAdd && (sampleId != sampleToAdd))
                    continue
                
                def sample_row = findSample(sampleId)
                if(!sample_row) {
                    sample_row = addSample(sampleId, peds, cohort, batch)
                }
                
                def variant_obs = db.firstRow("select * from variant_observation where sample_id = ${sample_row.id} and variant_id = ${variant_row.id} and batch_id = $batch;")
                if(!variant_obs) {
                    db.execute("""
                        insert into variant_observation (id,variant_id,sample_id, batch_id, qual,dosage, created) 
                                    values (NULL, $variant_row.id, 
                                                  ${sample_row.id}, 
                                                  ${batch},
                                                  ${v.sampleGenoType(sampleId)?.GQ?.toDouble()}, 
                                                  ${v.sampleDosage(sampleId)}, 
                                                  datetime('now'));
                    """)
                    ++countAdded
                }
            }
        }
        return countAdded
    }
    
    /**
     * Execute the given closure in teh scope of a transaction, and roll it back if
     * an exception occurs.
     * 
     * @param c     Closure to execute
     */
    void tx(Closure c) {
        db.execute("BEGIN TRANSACTION;")
        try {
            c()
            db.execute("COMMIT;")
        }
        catch(Exception e) {
            System.err.println "Database operation failed: $e"
            db.execute("ROLLBACK;")
            throw e
        }
    }
    
    void close() {
        if(this.db)
            this.db.close();
    }
    
    /**
     * Return a set of counts of times the given variant has been 
     * observed with in
     * 
     * a) all cohorts
     * b) all except the specified cohort
     * c) within the specified cohort
     */
    Map queryVariantCounts(Map options=[:], Variant variant, Allele allele, String sampleId, String cohort) {
        
        def params = [
            chr: variant.chr,
            pos: allele.start,
            alt: allele.alt,
            cohort: cohort
        ]

        String batchClause = ""
        if(options.batch) {
            batchClause = ' AND o.batch_id = :batch '
            params.batch = options.batch
        }
        
        // The total observations of the variant
        def variant_count = db.firstRow('''
            select count(distinct(o.sample_id))
            from variant_observation o, variant v 
            where o.variant_id = v.id and v.chr = :chr and v.pos = :pos and v.alt = :alt
        ''' + batchClause, params)[0]
        
        def excludeCohorts = ""
        if(options.excludeCohorts) 
            excludeCohorts = 'and s.cohort not in ("' + options.excludeCohorts.join('","') + '")'
            
       
        def variant_count_outside_target = 0
//        if(variant_count>1) { // Will always be at least 1 because it will be recorded for the sample we are processing now
            def outsideTargetSql = '''
            select count(distinct(s.id))
            from variant_observation o, variant v, sample s
            where o.variant_id = v.id 
                  and v.chr = :chr
                  and v.pos = :pos
                  and v.alt = :alt
                  and o.sample_id = s.id
                  and s.cohort <> :cohort
           ''' + excludeCohorts + batchClause
           
           println(outsideTargetSql)
            
            variant_count_outside_target = db.firstRow(outsideTargetSql, params)[0] 
            
        
            log.info "Variant $variant found $variant_count times globally, $variant_count_outside_target outside target $cohort / $options.excludeCohorts"
//        }
        
        int variant_count_within_target = db.firstRow('''
            select count(distinct(s.id))
            from variant_observation o, variant v, sample s
            where o.variant_id = v.id 
                  and v.chr = :chr
                  and v.pos = :pos
                  and v.alt = :alt
                  and o.sample_id = s.id
                  and s.cohort = :cohort
            ''' + batchClause, params)[0] 
         
        return [ total: variant_count, other_target: variant_count_outside_target, in_target: variant_count_within_target ]
    }
    
    /**
     * Simple test program - all it does is creates the database
     * 
     * @param args
     */
    static void main(String [] args) {
        println "Testing variant DB"
        def db = new VariantDB("test.db")
        println "Successfully connected to database"
    }
}
