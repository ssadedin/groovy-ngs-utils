package db

import groovy.sql.Sql
import groovy.util.logging.Log;

/**
 * The default schema, which uses SQLite syntax
 * 
 * @author simon.sadedin@mcri.edu.au
 */
@Log
class Schema {
    
    /**
     * Meta data about schema that this database is connected to
     */
    def schemaInfo = null
    
    
    static String CREATE_SAMPLE_SQL = 
           """
                create table sample (
                       id integer primary key asc, 
                       sample_id text NOT NULL,
                       father_id text,
                       mother_id text,
                       family_id text,
                       phenotype integer NOT NULL,
                       created datetime NOT NULL,
                       UNIQUE (sample_id) ON CONFLICT ROLLBACK
                );
           """
           
           
    static String VARIANT_OBSERVATION_TABLE = """ (
                       id integer primary key asc, 
                       variant_id integer references variant(id),
                       sample_id integer references sample(id),
                       batch_id text,
                       qual float,
                       dosage integer,  -- how many copies of the variant (1=het, 2=hom)
                       created datetime NOT NULL,
                       UNIQUE (variant_id,sample_id,batch_id) ON CONFLICT ROLLBACK
                )
    """

    static Map<Integer, Map<String,List<String>>>  schema = [
            1 : [ // Schema version 1
               """
                create table variant ( 
                        id integer primary key asc, 
                        chr text NOT NULL,
                        start integer,          -- the start of actual DNA change
                        end integer,            -- end end of actual DNA change
                        ref text NOT NULL,      -- reference sequence
                        alt text NOT NULL,      -- observed sequence
                        consequence text,       -- type of protein change (where applicable)
                        sift text,              -- sift impact
                        polyphen text,          -- polyphen impact
                        condel float,           -- condel sore 
                        max_freq float,         -- maximum frequency observed in population databases
                        dbsnp_id text           -- DBSNP ID (if known)
                );
            """,
            """
                CREATE INDEX variant_idx ON variant (chr,start,alt);
            """,
            """
                CREATE INDEX variant_cons_idx ON variant (consequence);
            """,
            """
                CREATE INDEX variant_condel_idx ON variant (condel);
            """,
            ' create table variant_observation ' + VARIANT_OBSERVATION_TABLE + ';',
           """
                CREATE INDEX variant_observation_idx ON variant_observation (variant_id);
           """,
           CREATE_SAMPLE_SQL,
           """
                CREATE INDEX sample_index ON sample (sample_id);
           """,
           """
                create table schema_meta_data (
                       id integer primary key asc, 
                       schema_version integer NOT NULL
                );    
           """,
           """
              insert into schema_meta_data values(1,1);
           """
          ],
      2: [
            """
              alter table sample add column cohort text 
              ;
            """,
            """
                CREATE INDEX sample_cohort_idx ON sample (cohort);
            """,
            """
              alter table sample add column batch text 
              ;
            """,
            """
                CREATE INDEX sample_batch_idx ON sample (batch);
            """
         ],
      3: [  
            """
              alter table variant add column pos integer 
            """
         ],
      4: [
          """
              alter table variant add column protein_change text;
          """
         ],
     5: [
         """
         create table sample_batch (
             id integer primary key asc,
             batch text,
             sample_id integer references sample(id),
             cohort text,
             created datetime NOT NULL
         );
         """,
         """
                CREATE INDEX sample_batch_sample_id_idx ON sample_batch (sample_id);
        """,
        """
                CREATE INDEX sample_batch_cohort_idx ON sample_batch (cohort);
        """
        ],
     6: [
            """
                    alter table variant add column gene text;
            """,

            """
                    CREATE INDEX gene_index ON variant (gene);
            """
        ]
    ]
    
    /**
     * Iterate over the schema versions and check if they exist
     */
    void checkSchema(Sql db) {
        
        
        // Get the current schema information
        try {
            schemaInfo = db.firstRow("select * from schema_meta_data")
        }
        catch(Exception e) {
            schemaInfo = [id:0, schema_version:0] // 0 meaning, "does not exist"
            
            boolean legacySchema = false
            try {
                db.firstRow("select count(*) from variant")
                
                // Variant table exists but not schema table
                legacySchema=true
            }
            catch(Exception e2) {
            }
            
            if(legacySchema) {
                upgradeLegacy(db)
                schemaInfo = db.firstRow("select * from schema_meta_data")
            }
        }
        
        List schemaVersions = (schema.keySet() as List).sort()
        schemaVersions.each { version ->
            if(schemaInfo.schema_version < version)  {
                upgradeSchema(db, version)
                schemaInfo.schema_version = version
            }
        } 
    }
    
    /**
     * This function upgrades a "legacy" database that was created as part of the
     * MGHA project to the schema used by this class. It's a one-way conversion, so 
     * use caution!
     * 
     * @param db
     */
    void upgradeLegacy(Sql db) {
        
        def statements = [
            """
                        alter table variant add column consequence text;
            """,
            """
                        alter table variant add column sift text;
            """,
            """
                        alter table variant add column polyphen text;
            """,
            """
                        alter table variant add column condel float;
            """,
            """
                        alter table variant add column max_freq float;
            """,
            """
                        update variant set max_freq = max(freq_1000g, freq_esp);
            """,
            """
                CREATE INDEX variant_max_freq_idx ON variant (max_freq);
            """,
            """
                ALTER TABLE sample RENAME TO tmp_sample;
            """,
            """
            create table sample (
                id integer primary key asc,
                sample_id text NOT NULL,
                father_id text,
                mother_id text,
                family_id text,
                phenotype integer NOT NULL,
                created datetime NOT NULL,
                UNIQUE (sample_id) ON CONFLICT ROLLBACK
            );
            """,
            """
              alter table sample add column cohort text 
              ;
            """,
            """
                CREATE INDEX sample_cohort_idx ON sample (cohort);
            """,
            """
              alter table sample add column batch text 
              ;
            """,
            """
                CREATE INDEX sample_batch_idx ON sample (batch);
            """,
          
            """
            INSERT INTO sample SELECT id, study_id, 0, 0, study_id, 0, created, cohort, batch   FROM tmp_sample;
            """,
            """
            DROP TABLE tmp_sample;
            """,
            """
                CREATE INDEX sample_index ON sample (sample_id);
            """, 
            """
                create table schema_meta_data (
                       id integer primary key asc, 
                       schema_version integer NOT NULL
                );    
            """,
            'create table variant_observation_tmp ' + VARIANT_OBSERVATION_TABLE + ';',
            """
              insert into variant_observation_tmp select id, variant_id, sample_id, NULL, qual,dosage,created from variant_observation;
            """,
            """
              drop table variant_observation;
            """,
            """
                alter table variant_observation_tmp rename to variant_observation;
            """,
            """
              insert into schema_meta_data values(1,4);
            """
        ]
                
        /*
        System.err.println("Legacy database upgrade required. The following statements will be executed:")
        
        for(s in statements) {
            System.err.println("="*80)
            System.err.println(s)
        }
        */
                
        for(s in statements) {
            System.err.println "Executing legacy upgrade statement: $s"
            db.execute(s)
        }
    }
    
    /**
     * Execute the given statements to upgrade the schema
     * 
     * @param statements
     */
    void upgradeSchema(Sql db, int toVersion) {
        log.info "Schema is out of date: upgrading to version $toVersion"
        db.execute("BEGIN TRANSACTION;")
        List<String> statements = schema[toVersion] 
        for(String s in statements) {
            log.info "Executing upgrade statement: $s"
            db.execute(s)
        }
        db.execute("update schema_meta_data set schema_version=$toVersion")
        log.info "Committing changes"
        db.execute("COMMIT;")
    }
}
