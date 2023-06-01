# Groovy NGS ![Tests](https://github.com/ssadedin/groovy-ngs-utils/actions/workflows/ci-build.yml/badge.svg)

A toolkit for working with genomic sequencing data in Groovy.

The JVM is an incredible platform for data analysis, offering high performance, extraordinary library and platform
support and rock solid industry support when it comes time to scale up and productionise your work.
Groovy NGS aims to unlock the power of the JVM for working with genomic sequencing data by enabling it 
to be used with the versatile and highly productive Groovy programming language.

Groovy NGS can be used at three levels:

 * Directly as pre-written tools on the command line
 * For writing simple scripts (bash-style) or interactive analysis in 
   [Jupyter Notebooks](https://github.com/ssadedin/beakerx)
 * As a library of classes for building full-scale applications
 
Under the hood, Groovy NGS is built on the widely used [HTSJDK](https://github.com/samtools/htsjdk). However
Groovy NGS makes it much easier to work with these libraries by adding idiomatic Groovy 
language constructs and filling in important commonly used missing features.
 
Examples of supported functionality are:

  * Reading, processing and filtering VCF files, including integration with common annotation sources such as VEP
  * Working with Genomic Ranges - full set of operation as well as higher level reading, processing and filtering 
  * Reading, processing and performing logical operations with pedigree (PED) files and family structures
  * Working with BAM/SAM/CRAM files (including, generating and working with Pileups)
  * A range of statistical operations including R-like data frames and linear modeling constructs
  * Many many more useful operations
 
For more information see:

- the [user guide](https://ssadedin.github.io/groovy-ngs-utils/userguide/) for general overview and examples
- the [reference documentation](http://ssadedin.github.io/groovy-ngs-utils/doc/overview-summary.html) for full details

## Build

Clone the repository:

```
    git clone --recursive git@github.com:ssadedin/groovy-ngs-utils.git
```

Run gradle:

```
    cd groovy-ngs-utils
    ./gradlew clean jar
```

Note: if behind a proxy, you can specify it like so:

```
./gradlew -Dhttp.proxyHost=<host> -Dhttp.proxyPort=<port> clean jar
```

## Install

GNGS doesn't need any installation itself. However, since it is using groovy,
you need Groovy installed. The groovy version is controlled in the 
`gradle.properties`. You should check the current version in there and ensure it
matches the version of groovy available in your environment.

If you don't have the right groovy version, you can install it easily
using [SDK man](http://sdkman.io/). With SDKMan, it is just:

```
sdk install groovy 3.0.10
```

Once groovy is available and the build instructions have worked, the `gngs` and `gngstool` scripts
should "just work":

```bash
# What type of variants are in my VCF?
./bin/gngs 'println(VCF.parse("some.vcf").countBy { it.type })'
```

You can get a GNGS enabled interactive Groovy Shell like this:

```
./bin/gngsh
groovy:000> new SAM("my.bam").basesAt("chr7", 117292917)
===> [A:5, total:5, C:0, T:0, D:0, G:0]
```

For the command line tools implemented by GNGS, you can run them using `gngstools`:

```
./bin/gngstool ExtractFASTQ -bam my.bam | bwa mem -p hg19.fasta - | samtools view -Sb - >  my.realigned.bam
```


## Command Line Scripting 


If you want to write command line scripts, the easiest way is using the `gngs` launcher
tool.  Here are some examples of simple command line scripts showing how GNGS can be used:

```bash
  # Get me a filtered VCF with QUAL>20 and DP > 5
  cat my.vcf | gngs 'VCF.filter { it.qual > 20 && it.info.DP.toInteger()>5 }' > filtered.vcf

  # Which read shave MAPQ = 0?
  cat my.bam | gngs 'SAM.eachRead { if(it.mappingQuality == 0) { println it.readName } }'
  
  # What's the median coverage of bedtools output?
  coverageBed -d  -abam test.bam -b test.bed | cut -f 6 | gngs 'Stats.read().median'

  # What are the bases piled up at 25870138?
  gngs 'println(new SAM("test.bam").basesAt("chr1", 25870138))'

  # What are the exons in the DVL1 gene?
  gngs 'println(RefGenes.download().getExons("DVL1").join(", "))'

  # What are the total number of bases overlapped by regions in my BED file?
  gngs 'println(new BED("./tests/data/small.overlaps.bed").load().reduce().size())'

```

These are only examples and barely scratch the surface of all the functions built into
groovy-ngs-utils. You can find documentation about the individual classes and methods
in the [API Documentation](http://ssadedin.github.io/groovy-ngs-utils/doc/index.html)

You can build the documentation with:

```
./gradlew groovydoc
```

## Programming 

If you want to write longer scripts or full programs using GNGS as a library, 
most of the classes live in the `gngs` package:
```
import gngs.*
```

Then nearly everything will be imported.

## How it Works
  
Everything is built upon Samtools, Picard Tools, BioJava and Apache commons-math. The jar file that 
is built bundles all the necessary libraries so that you can easily include them all with just one
classpath entry (or put it into your .groovy/lib).

Careful attention has been paid to make, wherever possible, operations operate on streaming data so that
memory is not a bottleneck in manipulating large data sets.

## Tools

GNGS also includes a bunch of pre-written tools which are mostly simple wrappers of 
the internal classes, but which expose them directly as command line tools in their 
own right. These are launched using the `gngstool` command: 

```
$ gngstool Sex test.vcf
================================================================================
Sex
================================================================================

MALE
 
```

Some useful tools include:

| Tool                   | Description                                                                            |
| ---------------------- | ---------------------------------------------------------------------------------------|
| VCFtoHTML              | Convert a VCF to HTML format, Diff VCFs                                                |
| Sex                    | Determine Sex from VCF, BAM or FASTQ file                                              |
| Table                  | Print out formated table from tab or comma separated data with many other features     |
| MeanCoverageEstimator  | A fast sampling based method to estimate mean coverage from a BAM file                 |
| Gaps                   | Identify regions of coverage below a threshold with statistics from BEDTools output    |

