groovy-ngs-utils
================

A collection of utilities for working with next generation (MPS) sequencing data in Groovy

This is a small collection of Groovy wrappers and interfaces that make it easy to perform 
scripting operations with Groovy to process NGS data.

The quality of this code is very alpha - it's simply what I've created as part of my PhD as I
do various data handling tasks. The kind of operations currently supported are:

  * Reading, processing and filtering VCF files
  * Reading, processing and filtering BED files
  * Working with SAM files (particularly, generating and working with Pileups)
  * Predicting Restriction Enzyme cut sites

Since these utilities are optimized for scripting, they live in the default Java package. This means you can 
easily write command line scripts, such as:

  cat my.vcf | groovy -e 'VCF.filter { it.qual > 20 && it.info.DP.toInteger()>5 }' > filtered.vcf

  cat my.bam | 'SAM.eachRead { if(it.mappingQuality == 0) { println it.readName } }'
  
These functions are all built upon Samtools, Picard Tools and BioJava. The jar file that is built bundles all
the necessary libraries so that you can easily include them all with just one classpath entry (or put it into
your .groovy/lib).

NOTE: To get the full functions you may need to copy sam.jar, picard-private-parts.jar, and tribble.jar 
(eg: from GATK) into the lib directory and build the library. Do this if you get ClassNotFound exceptions 
when running the code.
