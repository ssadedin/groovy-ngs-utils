

class Cli extends CliBuilder {
    public Cli() {
        this.writer = new PrintWriter(System.err)
    }
    
    OptionAccessor check(Object args, List required) {
        def opts = this.parse(args)
        def missing = required.grep { !opts[it] }
        if(missing) {
            System.err.println "The following options are required but not provided:\n"
            missing.each { 
                println "\t-$it"
            }
            this.usage()
            System.exit(1)
        }
    }
}
