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
import db.Schema;
import groovy.sql.Sql
import groovy.util.logging.Log;

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
        return findSample(sampleId)
    }
    
    /**
     * Find a variant in the database by its start position
     */
    def findVariant(Variant variant, Allele allele) {
        db.firstRow("select * from variant where chr=$variant.chr and start=$allele.start and alt=$allele.alt")
    }
    
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
           """select count(*) 
                       from variant_observation o,
                            variant v 
                       where o.variant_id = v.id
                         and v.chr = $v.chr
                         and v.start = $allele.start
                         and v.end = $allele.end
                         and v.alt = $allele.alt
        """)[0]
        
        int familyCount = db.firstRow(
           """select count(distinct(s.family_id))
                     from variant_observation o,
                          variant v, 
                          sample s
                     where o.variant_id = v.id
                         and v.chr = $v.chr
                         and v.start = $allele.start
                         and v.end = $allele.end
                         and v.alt = $allele.alt
                         and o.sample_id = s.id
        """)[0]
        
        return [samples: sampleCount, families: familyCount]
    }
    
    /**
     * Add the given variant to the database. If sampleToAdd is specified,
     * add it for only the given sample. Otherwise add it for all the samples that 
     * are genotyped to have the variant
     * 
     * @param annotations   optional annotations to draw from. These can be used when 
     *                      the annotations are not embedded in the VCF file. The only external
     *                      annotations supported at the moment are Annovar annotations.
     *                      The annotations can be any object having properties which 
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
                    
                    float onekgFreq = annotations["1000g2010nov_ALL"]?.isFloat() ? annotations["1000g2010nov_ALL"].toFloat() : 0.0f
                    float espFreq = annotations["ESP5400_ALL"]?.isFloat() ? annotations["ESP5400_ALL"].toFloat() : 0.0f
                     
                    float maxFreq = Math.max(onekgFreq, espFreq)
                    
                    db.execute("""insert into variant (id,chr,pos,start,end,ref,alt,consequence,protein_change,max_freq, dbsnp_id) 
                                   values (NULL, $v.chr, 
                                                 $v.pos, 
                                                 $allele.start, $allele.end, 
                                                 $v.ref, $allele.alt,
                                                 $annotations.ExonicFunc,
                                                 ${annotations.AAChange_RefSeq},
                                                 $maxFreq, 
                                                 ${annotations.dbSNP138})
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
                    sample_row = addSample(sampleId, peds, batch, cohort)
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
            e.printStackTrace()
            db.execute("ROLLBACK;")
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
    Map queryVariantCounts(Map options=[], Variant variant, Allele allele, String sampleId, String cohort) {
        
        // The total observations of the variant
        def variant_count = db.firstRow(
                """
        select count(*) 
            from variant_observation o, variant v 
            where o.variant_id = v.id and v.chr = $variant.chr and v.pos = $allele.start and v.alt = $allele.alt
        """)[0]
        
        def excludeCohorts = ""
        if(options.excludeCohorts) 
            excludeCohorts = /and s.cohort not in ("${options.excludeCohorts.join('","')}")/
        
        def variant_count_outside_target = 0
        if(variant_count>1) { // Will always be at least 1 because it will be recorded for the sample we are processing now
            variant_count_outside_target = db.firstRow("""
            select count(*) 
            from variant_observation o, variant v, sample s
            where o.variant_id = v.id 
                  and v.chr = $variant.chr 
                  and v.pos = $allele.start
                  and v.alt = $allele.alt
                  and o.sample_id = s.id
                  and s.cohort <> $cohort
                  """ + excludeCohorts
            )[0] 
        
            log.info "Variant $variant found $variant_count times globally, $variant_count_outside_target outside target $cohort / $options.excludeCohorts"
        }
        
        int variant_count_within_target = db.firstRow("""
            select count(*) 
            from variant_observation o, variant v, sample s
            where o.variant_id = v.id 
                  and v.chr = $variant.chr 
                  and v.pos = $allele.start
                  and v.alt = $allele.alt
                  and o.sample_id = s.id
                  and s.cohort = $target
            """)[0] 
         
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
