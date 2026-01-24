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

## Key Requirements

- The duplex trimmer must be able to process very large files in a reasonable time
- By far the most expensive step is the alignment part
- We can assume the system has many cores and large RAM
- To facilitate processing large files, a multithreaded architecture should be implemented
- It is critical that reads are written out in the same order that they are read in
- It is also important that the queuing overhead doesn't become the bottleneck, so
  make sure there is almost no overhead when processing a read that fails
  very quick tests to check if it is a candidate for trimming
- Based on these requirements, an appropriate multithreading architecture should be
  implemented
  
## Implementation Guidance

The duplex trimmer is to be implemented in Groovy as a GNGS tool, extending the `Toolbase` class.
It should use the `gngs.SAM` class to process a BAM file and identify duplex reads. The core
API to utilise is the `SAMRecord` class to retrieve read information.

## Implementation Plan

### Phase 1: Core Detection Logic
- [x] Create DuplexDetector class
  - [x] Method to extract soft clip information from CIGAR
  - [x] Method to calculate if soft clip length ≈ aligned length (within 10%)
  - [x] Method to extract aligned and soft clipped sequences
  - [x] Method to align first 50bp of soft clip vs last 50bp of aligned part using Align.local()
  - [x] Configurable alignment threshold (default 80% match)

### Phase 2: Tool Implementation
- [ ] Create DuplexTrimmer class extending ToolBase
  - [ ] CLI options: -i input, -o output, -t threshold, -a action, -s stats
  - [ ] Processing pipeline using SAM.eachRecord()
  - [ ] Write output using SAM.withWriter()
- [ ] Implement actions:
  - [ ] reject: Skip duplex reads entirely
  - [ ] trim: Remove soft clips, keep aligned portion
  - [ ] keep_best: Compare quality scores, keep higher quality portion

### Phase 3: Testing
- [x] Unit tests (DuplexDetectorTest.groovy)
  - [x] Normal read (no soft clips) → not duplex
  - [x] Soft clip too small → not duplex
  - [x] Soft clip wrong sequence → not duplex
  - [x] True duplex read → IS duplex
  - [x] Edge case: exactly 10% difference in lengths
  - [x] Hard clips should be ignored
- [ ] Integration tests (DuplexTrimmerTest.groovy)
  - [ ] Create test BAM files with mix of normal and duplex reads
  - [ ] Verify output BAM has correct number of reads
  - [ ] Verify duplex reads handled per selected action
  - [ ] Verify statistics are accurate

### Phase 4: Statistics and Reporting
- [ ] Track metrics:
  - [ ] Total reads processed
  - [ ] Duplex reads detected
  - [ ] Reads rejected/trimmed/kept
  - [ ] Average soft clip length
  - [ ] Average alignment score
- [ ] Output statistics in JSON or TSV format

## Design Notes
- Use Groovy 2.x compatible syntax: `[1,2,3] as byte[]` not `new byte[]{1,2,3}`
- Use 4 spaces for indentation, no tabs
- Default alignment threshold: 80% (40/50 bases match)
- Default length tolerance: 10%
- Check both leading and trailing soft clips
- Use `gngs.Align.local()` for sequence alignment
