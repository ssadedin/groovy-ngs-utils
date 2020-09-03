package gngs.tools
import org.codehaus.groovy.runtime.StackTraceUtils;

import db.VariantDB
import gngs.Cli
import gngs.Pedigrees
import gngs.ProgressCounter
import gngs.ToolBase
import gngs.VCF
import gngs.Variant
import groovy.util.logging.Log

/**
 * Utility to import variants to SQLite database rows
 * 
 * @author Simon Sadedin
 */
@Log
class VCFtoDB extends ToolBase {

    @Override
    public void run() {
        Pedigrees peds = opts.ped ? Pedigrees.parse(opts.ped) : null

        log.info "FOO"

        gngs.db.VariantDB db = new gngs.db.VariantDB(opts.db)
        def progress = new ProgressCounter(withRate:true)
        int exitCode = 0
        try {
            db.tx {
                VCF.parse(opts.vcf) { Variant v ->
                    // println "Adding variant $v to database"
                    db.add(opts.batch, peds, v)    
                    progress.count()
                }
            }
            progress.end()
        }
        catch(Exception e) {
            StackTraceUtils.sanitize(e).printStackTrace()
            exitCode = 1
        }
        finally {
            db.db.close()
        }

        System.exit(exitCode)
    }
    
    static void main(String[] args) {
        cli('VCFtoDB -vcf <VCF> -db <SQLite database file>', args) {
            vcf "VCF file to import to database", args:1, required:true
            ped "Pedigree file (PED format) for samples in VCF file", args:1
            db "SQLite database file to insert into", args:1, required:true
            batch "Name of batch to which sample belongs", args:1, required:true
        }
    }
}