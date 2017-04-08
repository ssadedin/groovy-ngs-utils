# Groovy NGS Utils

A collection of utilities for working with next generation (MPS) sequencing data in Groovy

This is a collection of Groovy wrappers and interfaces that make it easy to perform 
scripting operations with Groovy to process NGS data.

The kind of operations currently supported are:

  * Reading, processing and filtering VCF files, including integration with common annotation sources such as VEP, Annovar and SnpEFF
  * Reading, processing and filtering BED files or any source of genomic ranges
  * Reading, processing and performing logical operations with pedigree (PED) files and family structures
  * Working with SAM files (particularly, generating and working with Pileups)
  * Predicting Restriction Enzyme cut sites
  * A range of statistical operations including R-like data frames and linear modeling constructs

Since these utilities are optimized for scripting, they live in the default Java package. This means you can 
easily write command line scripts, such as:

```bash
  cat my.vcf | groovy -e 'VCF.filter { it.qual > 20 && it.info.DP.toInteger()>5 }' > filtered.vcf

  cat my.bam | groovy -e 'SAM.eachRead { if(it.mappingQuality == 0) { println it.readName } }'
  
  coverageBed -d  -abam test.bam -b test.bed | cut -f 6 | groovy -e 'Stats.read().median'
```
  
These functions are all built upon Samtools, Picard Tools, BioJava and Apache commons-math. The jar file that 
is built bundles all the necessary libraries so that you can easily include them all with just one
classpath entry (or put it into your .groovy/lib).

Careful attention has been paid to make, wherever possible, operations operate on streaming data so that
memory is not a bottleneck in manipulating large data sets.

## Building

Clone the repository:

    git clone git@github.com:ssadedin/groovy-ngs-utils.git

Run gradle:

    cd groovy-ngs-utils
    ./gradlew jar


## Installation

The easiest way is to put the file into your ~/.groovy/lib folder:

```bash
    mkdir -p ~/.groovy/lib
    cp build/libs/groovy-ngs-utils.jar ~/.groovy/lib
```

This way you can use commands without modifying your classpath at all.

However: since this places groovy-ngs-utils.jar into your default groovy path
it is probably not a great idea since it may introduce versions of classes to your default Groovy
classpath that conflict with classes other applications using Groovy depend on.
It's an easy way to experiment with the library but you are better off
specifying it explicitly when you need it:

```bash
# What type of variants are in my VCF?
groovy -cp ~/groovy-ngs-utils/build/libs/groovy-ngs-utils.jar \
       -e 'println(VCF.parse("some.vcf").countBy { it.type })'
```
