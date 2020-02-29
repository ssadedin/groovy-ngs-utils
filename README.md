# Groovy NGS Utils [![travis](https://travis-ci.org/ssadedin/groovy-ngs-utils.svg?branch=master)](https://travis-ci.org/ssadedin/groovy-ngs-utils)

A library for working with next generation (MPS) sequencing data in Groovy.

The kind of operations supported are:

  * Reading, processing and filtering VCF files, including integration with common annotation sources such as VEP, Annovar and SnpEFF
  * Reading, processing and filtering BED files or any source of genomic ranges
  * Reading, processing and performing logical operations with pedigree (PED) files and family structures
  * Working with BAM/SAM/CRAM files (particularly, generating and working with Pileups)
  * Predicting Restriction Enzyme cut sites
  * A range of statistical operations including R-like data frames and linear modeling constructs

The code in this library is usable in three different ways:

 * Directly as tools on the command line
 * For writing simple scripts (bash-style)
 * As a library of classes for building full-scale applications


## Build

Clone the repository:

```
    git clone git@github.com:ssadedin/groovy-ngs-utils.git
    git submodule update --init
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

GNGS doesn't actually need any installation. However, since it is using groovy,
you need Groovy 2.4.x installed. Although other 2.4 versions should
work, I suggest using 2.4.10 since that version is what GNGS is tested
with currently. You can install it easily without any administrative
using [SDK man](http://sdkman.io/). With SDKMan, it is just:

```
sdk install groovy 2.4.10
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

