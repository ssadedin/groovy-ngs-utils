# Groovy NGS Utils

A collection of utilities for working with next generation (MPS) sequencing data in Groovy

This is a collection of Groovy wrappers and interfaces that make it easy to perform 
scripting operations with Groovy to process NGS data.

The quality of this code is beta level - it's simply what I've created as part of my PhD as I
do various data handling tasks. There are a decent number of unit tests, but they probably do not
achieve more than 50% code coverage. 

The kind of operations currently supported are:

  * Reading, processing and filtering VCF files, including integration with common annotation sources such as VEP, Annovar and SnpEFF
  * Reading, processing and filtering BED files or any source of genomic ranges
  * Reading, processing and performing logical operations with pedigree (PED) files and family structures
  * Working with SAM files (particularly, generating and working with Pileups)
  * Predicting Restriction Enzyme cut sites
  * A range of statistical operations including R-like data frames and linear modeling constructs

Since these utilities are optimized for scripting, they live in the default Java package. This means you can 
easily write command line scripts, such as:

  cat my.vcf | groovy -e 'VCF.filter { it.qual > 20 && it.info.DP.toInteger()>5 }' > filtered.vcf

  cat my.bam | groovy -e 'SAM.eachRead { if(it.mappingQuality == 0) { println it.readName } }'
  
  coverageBed -d  -abam test.bam -b test.bed | cut -f 6 | 'Stats.read().median'
  
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
    ./gradlew


## Installation

The easiest way is to put the file into your ~/.groovy/lib folder:

    mkdir -p ~/.groovy/lib
    cp build/libs/groovy-ngs-utils.jar ~/.groovy/lib

This way you can use commands without modifying your classpath at all.

*Note:* In rare circumstances placing a large Groovy library like this into
your default .groovy/lib path *could* cause some incompatibilities with other
libraries or tools if they were built with a different version of Groovy. This
sort of problem is rare, but worth keeping in mind if you find other 
Groovy based tools give strange errors.
