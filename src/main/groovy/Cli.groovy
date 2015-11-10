import org.apache.commons.cli.Option

class Cli extends CliBuilder {
    
    // To avoid script code having to do the above import, make an alias for it
    static int UNLIMITED = Option.UNLIMITED_VALUES
    
    public Cli() {
        this.writer = new PrintWriter(System.err)
    }
    
	void banner(String title) {
        System.err.println "="*100
		System.err.println " $title ".center(100)
        System.err.println "="*100
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
