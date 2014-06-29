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
                        consequence text,       -- protein change (where applicable)
                        sift text,              -- sift impact
                        polyphen text,          -- polyphen impact
                        condel float,           -- condel sore 
                        max_freq float,         -- maximum frequency observed in population databases
                        dbsnp_id text           -- DBSNP ID (if known)
                );
            """,
            """
                CREATE INDEX variant_idx ON variant (chr,start,alt);
                CREATE INDEX variant_cons_idx ON variant (consequence);
                CREATE INDEX variant_condel_idx ON variant (condel);
            """,
            """
                create table variant_observation (
                       id integer primary key asc, 
                       variant_id integer references variant(id),
                       sample_id integer references sample(id),
                       batch_id text,
                       qual float,
                       dosage integer,  -- how many copies of the variant (1=het, 2=hom)
                       created datetime NOT NULL,
                       UNIQUE (variant_id, sample_id) ON CONFLICT ROLLBACK
                );
           """,
           """
                CREATE INDEX variant_observation_idx ON variant_observation (variant_id);
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
     * Execute the given statements to upgrade the schema
     * 
     * @param statements
     */
    void upgradeSchema(Sql db, int toVersion) {
        log.info "Schema is out of date: upgrading to version $toVersion"
        db.execute("BEGIN TRANSACTION;")
        List<String> statements = schema[toVersion] 
        for(String s in statements) {
            db.execute(s)
        }
        db.execute("COMMIT;")
    }
}
