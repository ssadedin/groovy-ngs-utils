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
import java.sql.Timestamp;

import com.xlson.groovycsv.PropertyMapper;

import db.Schema;
import groovy.sql.Sql
import groovy.util.logging.Log;

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
class VariantDB implements Closeable {
    
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
                                ${family?.id}, 
                                ${subject?.phenoTypes?.getAt(0)?:0}, 
                                ${cohort},
                                ${batch},
                                datetime('now'));
        """)
        
        def result =  findSample(sampleId)
        addSampleBatch(result.id,batch,cohort)
        return result
    }
    
    void addSampleBatch(Long sampleId, String batch, String cohort) {
        db.execute("""
            insert into sample_batch (id, sample_id, batch, cohort, created) 
                 select NULL, 
                        $sampleId,
                        $batch,
                        $cohort,
                        datetime('now')
                 where not exists (select 1 from sample_batch where sample_id = $sampleId and batch = $batch);
        """)
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
     * Two values are returned: 
     *   <li>sampleCount - the number of unique samples in which the variant was observed
     *   <li>familyCount - the number of unique families in which the variant was observed
     *   
     * @param batch If provided, variants will only be provided if they are observed in
     *              this batch OR any earlier batch. The intent is to support reproducibility
     *              so that new samples can be added to the database without changing the
     *              variant counts if queried for old samples.
     */
    Map countObservations(Variant v, Allele allele, String batch=null) {
        
        def dateString = queryBatchDateTime(batch)
        
        int sampleCount = db.firstRow(
           """select count(distinct(o.sample_id))
                       from variant_observation o,
                            variant v 
                       where o.variant_id = v.id
                         and v.chr = $v.chr
                         and v.start = $allele.start
                         and v.end = $allele.end
                         and v.alt = $allele.alt
                         and o.created <= datetime($dateString)
        """)[0]
        
        int familyCount = db.firstRow(
           """select count(distinct(coalesce(s.family_id,s.sample_id)))
                     from variant_observation o,
                          variant v, 
                          sample s
                     where o.variant_id = v.id
                         and v.chr = $v.chr
                         and v.start = $allele.start
                         and v.end = $allele.end
                         and v.alt = $allele.alt
                         and o.sample_id = s.id
                         and o.created <= datetime($dateString)
        """)[0]
        
        return [samples: sampleCount, families: familyCount]
    }

    /**
     * Find the time of the most recent variant observation added for the given
     * batch, formatted correctly as a string relative to GMT for querying via SQLite. 
     * If batch is null, return the current time plus 1 second.
     * 
     * @param batch id of batch to query
     * @return a String containing the date in SQLite format
     */
    private String queryBatchDateTime(String batch) {
        def dateString = null
        if(batch) {
            dateString = db.firstRow("""
                select max(vo.created) 
                       from variant_observation vo 
                       where vo.batch_id = $batch
            """)[0]
        }
        else {
            Timestamp countBefore = new Timestamp(System.currentTimeMillis()+1000)
            dateString = countBefore.format("yyyy-MM-dd HH:mm:ss", TimeZone.getTimeZone("UTC"))
        }
        return dateString
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
     * Map of string to database rows of cached sample information
     */
    Map<String, Object> cachedSampleInfo = [:]
    
    /**
     * Returns cached rows from the sample table for each sample. If rows do not exist yet for
     * the sample, adds rows to the table. If rows do exist, returns the existing rows.
     * Also adds batch information.
     * 
     * @param addSamples
     * @param batch
     * @param cohort
     * @return a map keyed on sample id, with values being the row from the sample 
     *         table for each sample
     */
    Map addCachedSampleInfo(List<String> addSamples, String batch, String cohort, Pedigrees peds=null) {
        addSamples.collectEntries { sampleId ->
            def sample_row = cachedSampleInfo[sampleId]
            if(!sample_row) {
                sample_row = findSample(sampleId)
                if(!sample_row) {
                    sample_row = addSample(sampleId, peds, cohort, batch)
                }
                else {
                    if(batch != sample_row.batch) {
                        addSampleBatch(sample_row.id, batch, cohort)
                    }
                }
            }
            [sampleId, sample_row]
        }
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
        
        List addSamples = sampleToAdd ? [sampleToAdd] : v.header.samples
        Map sampleInfo = addCachedSampleInfo(addSamples, batch, cohort, peds)
        
        int countAdded = 0
        def alleles = alleleToAdd ? [alleleToAdd] : v.alleles
        for(Allele allele in v.alleles) {
            def variant_row = findVariant(v, allele)
            if(!variant_row) {
                
                if(annotations == null && v.header.hasInfo("CSQ")) {
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
                    
                    if(annotations == null)
                        annotations = [:]
                    
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
                                                 ${getAnnotation("DBSNP", annotations)},
										         $gene
                                                )
                               """)
                }
            }
            variant_row = findVariant(v,allele)
            
            // For every sample carrying the allele, add an observation
            for(String sampleId in addSamples.grep { v.sampleDosage(it) }) {
                
                def sample_row = sampleInfo[sampleId]
                
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
     * <li> all cohorts
     * <li> all except the specified cohort
     * <li> within the specified cohort
     * 
     * Additional cohorts can be excluded by setting a comma separated list as the
     * value 'excludeCohorts' in options.
     * <p>
     * The returned map has 3 attributes:
     * <li>total - total number of smaples observed to have the variant, regardless of cohort
     * <li>in_target - total number of samples observed to have the variant within specified cohort
     * <li>other_target - total number of samples observed to have the variant outside of specified cohort
     */
    Map queryVariantCounts(Map options=[:], Variant variant, Allele allele, String sampleId) {
        
        String dateString = queryBatchDateTime(options.batch)
        
        // The total observations of the variant
        def variant_count = db.firstRow(
            """
            select count(distinct(o.sample_id)) 
                from variant_observation o, variant v 
                where o.variant_id = v.id 
                  and v.chr = $variant.chr 
                  and v.pos = $allele.start 
                  and v.alt = $allele.alt
                  and o.created <= datetime($dateString)
            """)[0]
            
        def excludeCohortsInClause = ""
        
        List excludedCohorts = []
        if(options.excludeCohorts) 
            excludedCohorts.addAll(options.excludeCohorts)
        
        if(sampleId != null) {
            def sample_batch = db.firstRow("select coalesce(sb.cohort,s.cohort) as cohort from sample s, sample_batch sb where sb.sample_id = s.id and s.sample_id = $sampleId")
            def sampleCohorts = sample_batch?.cohort
            if(sampleCohorts) {
                log.info "Excluding cohorts $sampleCohorts based on cohorts of $sampleId"
                excludedCohorts.addAll(sampleCohorts)
            }
        }
        
        if(excludedCohorts)
            excludeCohortsInClause = 
              'and not exists (select 1 from sample_batch sb, sample s2 where s2.id = s.id and sb.cohort in ("' + excludedCohorts.join('","') + '"))'

          println "Running query: " +  """
            select count(distinct(s.id))
            from variant_observation o, variant v, sample s
            where o.variant_id = v.id 
                  and v.chr = $variant.chr 
                  and v.pos = $allele.start
                  and v.alt = $allele.alt
                  and o.sample_id = s.id
                  and o.created <= datetime($dateString)
                  """ + excludeCohortsInClause

              
                      
        def variant_count_outside_target = 0
        Utils.time("query_count_outside_target") {
            variant_count_outside_target = db.firstRow("""
            select count(distinct(s.id))
            from variant_observation o, variant v, sample s
            where o.variant_id = v.id 
                  and v.chr = $variant.chr 
                  and v.pos = $allele.start
                  and v.alt = $allele.alt
                  and o.sample_id = s.id
                  and o.created <= datetime($dateString)
                  """ + excludeCohortsInClause
            )[0] 
        }
        
        log.info "Variant $variant found in $sampleId $variant_count times globally, $variant_count_outside_target outside $excludedCohorts / $options.excludeCohorts"
        
        log.info "Running query: " + """
            select count(distinct(vs.id))
            from variant_observation o, variant v, sample qs, sample vs, sample_batch qsb, sample_batch vsb
            where o.variant_id = v.id 
                  and v.chr = $variant.chr 
                  and v.pos = $allele.start
                  and v.alt = $allele.alt
                  and o.sample_id = vs.id
                  and vsb.sample_id = vs.id
                  and qs.sample_id = $sampleId
                  and qsb.sample_id = qs.id
                  and qs.cohort = vs.cohort
                  and o.created <= datetime($dateString)
        """
        
        int variant_count_within_target = Utils.time("query_count_within_target") {
             db.firstRow("""
                select count(distinct(vs.id))
                from variant_observation o, variant v, sample qs, sample vs, sample_batch qsb, sample_batch vsb
                where o.variant_id = v.id 
                      and v.chr = $variant.chr 
                      and v.pos = $allele.start
                      and v.alt = $allele.alt
                      and o.sample_id = vs.id
                      and vsb.sample_id = vs.id
                      and qs.sample_id = $sampleId
                      and qsb.sample_id = qs.id
                      and qs.cohort = vs.cohort
                      and o.created <= datetime($dateString)
                """)[0] 
        }
             
        return [ 
            total: variant_count, 
            other_target: variant_count_outside_target, 
            in_target: variant_count_within_target 
        ]
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
