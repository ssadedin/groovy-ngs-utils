# Groovy NGS

A toolkit for working with genomic sequencing data in Groovy.

The JVM is a powerful platform for data analysis, offering high performance, strong library and platform
support with excellent deployment options when it comes time to scale up and productionise your work.
However Java itself is often a cumbersome language to work with. Groovy NGS aims to give the best 
of both worlds - the power of the JVM, combined with the productivity of a flexible and user
friendly scripting language.
 
Under the hood, Groovy NGS is built on the widely used [HTSJDK](https://github.com/samtools/htsjdk). However
Groovy NGS makes it much easier to work with these libraries by adding idiomatic Groovy 
language constructs and filling in important commonly used missing features.

Groovy NGS can be used at three levels:

 * Directly as pre-written tools on the command line
 * For writing simple scripts (bash-style) or interactive analysis in 
   [Jupyter Notebooks](https://github.com/ssadedin/beakerx)
 * As a library of classes for building full-scale applications
 
Examples of supported functionality are:

  * Reading, processing and filtering VCF files, including integration with common annotation sources such as VEP
  * Working with Genomic Ranges - full set of operation as well as higher level reading, processing and filtering 
  * Reading, processing and performing logical operations with pedigree (PED) files and family structures
  * Working with BAM/SAM/CRAM files (including, generating and working with Pileups)
  * A range of statistical operations including R-like data frames and linear modeling constructs
  * Many many more useful operations
 
