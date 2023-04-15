# As a Dependency In Projects

The default way to use Groovy NGS is to build a "fat" jar and include it in your project.
This jar file can be created by cloning the repository and executing the gradle command
to build the jar:

```
git clone --recursive https://github.com/ssadedin/groovy-ngs-utils.git
./gradlew jar
```

The jar file is created in `build/libs/groovy-ngs-utils.jar` and includes all necessary dependencies
*except* for Groovy itself. Therefore, you should include the [groovy-all](https://mvnrepository.com/artifact/org.apache.groovy/groovy-all)
dependency from Maven for the appropriate Groovy version in your project as well.

Note: in the future it is expected that Groovy NGS will move to support Maven 
repository based dependency resolution so that it can be more easily included in
other projects.

# On the Command Line

Another way to use Groovy NGS is to use the `gngs` command to execute expressions directly
on the command line. For example, after building (as above), to filter a VCF to only chromosome 21:

```bash
cat test.vcf | ./bin/gngs 'VCF.filter { it.chr == "chr21" }' > test.chr21.vcf
```

# Interactively via a Groovy Shell

You can run Groovy NGS as a [REPL](https://en.wikipedia.org/wiki/Read%E2%80%93eval%E2%80%93print_loop) to run commands
interactively and see the results using the standard Groovy Shell. To facilitate this, there is a script in the `bin` 
folder to launch the Groovy Shell with Groovy NGS added to the classpath and the gngs package automatically imported. This
can be a very useful way to experiment and learn about how Groovy NGS and Groovy itself works.

Here is how an example session looks:

```bash
./bin/gngsh 
Groovy Shell (3.0.10, JVM: 11.0.18)
Type ':help' or ':h' for help.
------------------------------------------------------------------------------------------------------
groovy:000> import gngs.*;
===> gngs.*
groovy:000> vcf = VCF.parse('src/test/data/giab1.tiny.trio.vcf')
===> VCF file src/test/data/giab1.tiny.trio.vcf for NA12877,NA12878,NA12879
groovy:000> vcf.size()
===> 24
groovy:000> vcf.count { it.het }
===> 20
```


# In a Jupyter Notebook

A final way to use Groovy NGS is via Jupyter. You can use any available Groovy kernel, however a 
particularly useful one is [BeakerX](https://github.com/twosigma/beakerx) which provides Groovy 
kernels out of the box along with a range of Jupyter widgets as enhancements. Together with Groovy NGS.
this makes a very effective data analysis platform for genomic data.



