package gngs.tools

import java.util.regex.Pattern

import gngs.Region
import gngs.SAM
import gngs.ToolBase
import groovy.json.JsonBuilder
import groovy.transform.CompileStatic
import groovy.transform.ToString
import groovy.util.logging.Log
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SAMTag

/**
 * Extracts a summary of the secondary alignments for reads in a specified region.
 * 
 * Outputs a list of JSON elements of the form:
 * 
 * <pre>
 *     {
 *      "readid": "1d6734cc-babd-42a1-9197-67b40a5ddcf8",
 *      "chrom": "chr21",
 *      "pos": 13836134,
 *      "strand": "+",
 *      "qual": 60,
 *      "offset": 28,
 *      "rlen": 22324,
 *      "qlen": 22283
 *  }
 * </pre>
 * Used for plotting structural rearrangemnts
 */
@Log
class ExtractSecondaryAlignmentSegments extends ToolBase {

    
    /**
     * Scan reads from given BAM / region to extract those with secondary alignments
     * and break the secondary alignments into separate SASegment objects
     * 
     * @param bam
     * @param region
     * @return List of tuples of form (readName, List of SASegments)
     */
    @CompileStatic
    List<List> scanReads(SAM bam, Region region) {
        
        return (List<List>)bam.withIterator(region) {  i ->
            i.grep { SAMRecord r ->
                !r.readUnmappedFlag && r.hasAttribute(SAMTag.SA)
            }.collect { o ->
                
                SAMRecord r  = (SAMRecord) o
                String strand = r.readNegativeStrandFlag ? "-" : "+"
                SASegment primaryAlignment = 
                    new SASegment(chr: r.referenceName, pos: r.alignmentStart, strand: strand, cigar: r.cigarString, qual: r.mappingQuality)

                List<SASegment> segs = 
                    [primaryAlignment] +

                    r.getStringAttribute(SAMTag.SA)
                     .split(";")
                     .grep { it }
                     *.split(",")
                     .collect { String[] parts ->
                        new SASegment(
                            chr: parts[0],
                            pos: parts[1].toInteger(),
                            strand: parts[2],
                            cigar: parts[3],
                            qual: parts[4].toInteger()
                        )
                    }
               return [r.readName, segs]
            }
        }

    }

    @Override
    public void run() {
        
        Region locus = new Region(opts.locus)
        
        log.info "Scanning locus $locus in $opts.bam"

        List<List> reads = scanReads(new SAM(opts.bam), locus)
        
        log.info "Finsihed scan: found ${reads.size()} reads"
        
       
        List<Map> results = reads.collectMany { List item ->
            
           String name = item[0]

            return item[1].collect { SASegment seg ->

                Map cigarSplit = seg.splitCigar()
               
                return [
                    "readid": name,
                    "chrom": seg.chr,
                    "pos": seg.pos,
                    "strand": seg.strand,
                    "qual": seg.qual,
                    * : cigarSplit
                ]
            }
            
        }
        println(new JsonBuilder(results).toPrettyString())
    }
    
    static void main(String [] args) {
        cli('Extract alignment segments', args) {
            bam 'BAM File', args:1, required: true
            locus 'Locus to search', args:1, required: true, type: String
        }
    }
}

/**
 * A supplementary alignment segment
 */
@ToString(includeNames=true)
class SASegment {
    String chr
    int pos
    String strand
    String cigar
    int qual
    
    final static Pattern cigarSplitter = ~"([0-9]+)([MIDSH])"
    
    /**
     * Compute the total query length, mapped length and offset (based on soft clips)
     * of the CIGAR string of this segment
     * 
     * @return  Map with rlen, qlen and offset attributes
     */
    Map<String, Integer> splitCigar() {
        
        List<String> cig =
                cigarSplitter.matcher(this.cigar).collect { [it[1].toInteger(),it[2]] }

        if(this.strand == "-")
            cig = cig.reverse()

        return [
            // sum length of all the CIGAR elements that represent reference sequence span
            rlen : cig.grep { it[1] == 'M' || it[1] == 'D' }.sum { it[0] }?:0,

            // sum length of all the CIGAR elements that represent query sequence span
            qlen  : cig.grep { it[1] == 'M' || it[1] == 'I' }.sum { it[0] }?:0,

            // length of any leading clipped sequence
            offset : cig[0][1] in ['S','H'] ? cig[0][0] : 0
        ]
     }
}
