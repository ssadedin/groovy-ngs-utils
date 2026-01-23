# Duplex Trimming

The Duplex Trimmer is a tool that detects reads that were accidentally read in duplex form
but were processed as simplex reads by an Oxford Nanopore instrument. When a read is read
in duplex form, it is read twice because the the first strand is read and then the
second strand of the same molecule is captured by the pore and read as a continuous
fragment. The duplex read has the following characteristics:

- the read is comprised of two key parts:
  - an aligned part
  - a soft clipped part
- the aligned part is approximately the same length as the soft clipped part
- the soft clipped part has approximately the same sequence as the aligned part
  but reverse complemented
  
 To identify the duplex read, one would:
 
- check if the total length of the soft clipped sequence was within 10% of the length of the total
  number of bases in the read
- extract the sequence of aligned bases and the soft clipped bases
- try to align the first 50 bases of the of the soft clipped part against
  the last 50 bases of the non-soft clipped part
  - if alignment exceeds a specified threshold, flag the read as a duplex read

When a read is flagged as a duplex read, there could be several options. One option is simply
to reject that read altogether. A second option would be to split off the soft clips so that
only the non-soft clips are retained.  A third option would be to split the read into two separate
reads (this is not advisable because it would align identically and behave as a duplicate).
Finally, it could be that we decide which of the two has higher quality overall 
and retain only the portion with the highest quality.

## Implementation Guidance

The duplex trimmer is to be implemented in Groovy as a GNGS tool, extending the `Toolbase` class.
It should use the `gngs.SAM` class to process a BAM file and identify duplex reads. The core
API to utilise is the `SAMRecord` class to retrieve read information.