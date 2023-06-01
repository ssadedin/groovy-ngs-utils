# Groovy NGS

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
 
