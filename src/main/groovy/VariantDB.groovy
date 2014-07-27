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
    def addSample(String sampleId, Pedigrees peds) {
        Pedigree family = peds.subjects[sampleId]
        Subject subject = family?.individuals.find { it.id == sampleId }
        if(!subject) 
            throw new IllegalStateException("Sample $sampleId could not be located in the pedigree file. Known samples are ${peds.subjects*.key}")
            
        sampleId = trimSampleId(sampleId)
                
        db.execute("""
            insert into sample (id, sample_id, father_id, mother_id, family_id, phenotype, created) 
                        values (NULL, 
                                $sampleId, 
                                ${family?.motherOf(sampleId)?.id}, 
                                ${family?.fatherOf(sampleId)?.id}, 
                                ${family?.id}, 
                                ${subject?.phenoTypes?.getAt(0)}, 
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
                         and o.sample_id = s.id
        """)[0]
        
        return [samples: sampleCount, families: familyCount]
    }
    
    /**
     * Add the given variant to the database
     */
    void add(String batch, Pedigrees peds, Variant v) {
        
        for(Allele allele in v.alleles) {
            def variant_row = findVariant(v, allele)
            if(!variant_row) {
                
                def vep = v.getVepInfo()[allele.index]
                def cons = v.getConsequence(allele.index)
                if(cons != null && cons.indexOf("(")>=0) 
                    cons = cons.substring(0, rawCons.indexOf("("))
                db.execute("""
                    insert into variant (id,chr,start,end,ref,alt, sift, polyphen, condel, consequence, max_freq,dbsnp_id) 
                                values (NULL, $v.chr, $allele.start, $allele.end, ${v.ref}, $allele.alt, ${vep?.SIFT}, ${vep?.PolyPhen}, ${vep?.Condel}, 
                                       ${cons}, $v.maxVepMaf, $v.id);
                """)
            }
            variant_row = findVariant(v,allele)
            
            // For every sample carrying the allele, add an observation
            for(String sampleId in v.header.samples.grep { v.sampleDosage(it) }) {
                
                def sample_row = findSample(sampleId)
                if(!sample_row) {
                    sample_row = addSample(sampleId, peds)
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
                }
            }
        }
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
