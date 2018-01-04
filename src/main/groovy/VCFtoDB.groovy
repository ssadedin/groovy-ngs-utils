import org.codehaus.groovy.runtime.StackTraceUtils;

import gngs.Cli
import gngs.Pedigrees
import gngs.ProgressCounter
import gngs.VCF
import gngs.Variant

Cli cli = new Cli()
cli.with {
    vcf "VCF file to import to database", args:1, required:true
    ped "Pedigree file (PED format) for samples in VCF file", args:1, required:true
    db "SQLite database file to insert into", args:1, required:true
    batch "Name of batch to which sample belongs", args:1, required:true
}

def opts = cli.parse(args)
if(!opts)
    System.exit(1)

Pedigrees peds = Pedigrees.parse(opts.ped)

VariantDB db = new VariantDB(opts.db)
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