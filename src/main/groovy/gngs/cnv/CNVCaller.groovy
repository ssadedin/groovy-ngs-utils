package gngs.cnv

import java.io.Writer

import org.apache.commons.math3.distribution.NormalDistribution
import gngs.GRange
import gngs.RangeIndex
import gngs.Regions
import gngs.RegulatingActor
import gngs.Utils
import graxxia.ThresholdCondition
import graxxia.ThresholdRange
import graxxia.Thresholder
import groovy.transform.CompileStatic
import groovy.util.logging.Log

/**
 * Processes {@link Contrast} objects and detects those that show a significant
 * difference in coverage between the detection region and the non-detection regions.
 * <p>
 * Significant regions are expanded to their corresponding window size and accumulated in the
 * {@link #cnvs} object.
 *  
 * @author Simon Sadedin
 */
@Log
@CompileStatic
class CNVCaller extends RegulatingActor<Contrast> implements ThresholdCondition {
    
    private List<Integer> widths
    
    Writer errorOut
    Writer countsOut
    Writer likelihoodsOut
    
    boolean headerWritten = false

    List<Thresholder> thresholders  = null
    
    /**
     * The CNVs called, stratified by size
     */
    Map<Integer, RangeIndex> cnvs = [:]
    
    /**
     * Create a CNV caller actor
     * 
     * @param chr               the chromosome that this caller will process
     * @param errorOut          
     * @param countsOut
     * @param likelihoodsOut
     */
    public CNVCaller(String chr, Writer errorOut, Writer countsOut, Writer likelihoodsOut) {
        super(10000,50000)
        this.errorOut = errorOut;
        this.countsOut = countsOut;
        this.likelihoodsOut = likelihoodsOut;
    }
    

    @Override
    public void process(Contrast contrast) {
        if(!headerWritten) {
            writeHeaders(contrast)
       }
        
        updateLikelihoods(contrast)
        
        writeErrors(contrast)
        
        writeCounts(contrast)
    }

    private void writeHeaders(Contrast contrast) {
        
         widths = (List<Integer>)contrast.layers.counts*.widths.flatten()
        
        // note: subtracting 1 from width is a hack to restore the original width because it gets expanded
        // by 1 through the range operations
        if(errorOut != null) {
            String errorHeader = (
                    ['pos'] + (List)contrast.layers.counts*.collect { LayerCount lc -> 'w'+(lc.widths[1]-1) }*.getAt(0)
                    ).join('\t')
            errorOut.write(errorHeader)
            errorOut.write('\n')
            log.info "Wrote error header: " + errorHeader
        }

        if(countsOut != null) {
            String countHeader = (
            ['pos', *contrast.layers.counts*.collect { LayerCount lc ->
                    int w = lc.widths[1] -1
                    [ 'L'+ w, 'M'+w, 'R' + w ]
                }*.getAt(0).flatten()
            ]).join('\t')
            log.info "Wrote counts header: " + countHeader
    
            countsOut.write(countHeader)
            countsOut.write('\n')
        }

        likelihoodsOut.write(['pos', *(contrast.layers.counts*.collect { LayerCount lc ->
                'w' + (lc.widths[1] -1)
            }*.getAt(0))].join('\t') + '\n')
        headerWritten = true
    }
    
    private final static char TAB_CHAR = '\t' as char

    private void writeCounts(final Contrast contrast) {
        
        if(countsOut.is(null))
            return
        
        final StringBuilder line = new StringBuilder()
        line.append(contrast.position)
        final LayerCount [] counts = contrast.layers.counts
        for(int i=0; i<counts.length; ++i) {
            final int[] count = counts[i].counts
            for(int j=0; j<count.length; ++j) {
                line.append(TAB_CHAR)
                line.append(String.valueOf(count[j]))
            }
        }
        line.append('\n')
        countsOut.write(line.toString())
    }

    private void writeErrors(Contrast contrast) {
        if(errorOut.is(null))
            return
        errorOut.write(((List)[contrast.position] + contrast.squaredErrors).join('\t') + '\n')
    }

    private void updateLikelihoods(Contrast contrast) {
        final List<Double> deletionLRs = (List<Double>)contrast.deletionLikelihoods

        updateThresholders(contrast, deletionLRs)
        likelihoodsOut.write(String.valueOf(contrast.position))
        
        final int numLRs = deletionLRs.size()
        for(int i=0; i<numLRs; ++i) {
            likelihoodsOut.write('\t')   
            likelihoodsOut.write(String.valueOf(Math.round(deletionLRs.getAt(i))))
        }
        likelihoodsOut.write('\n')
    }
    
    boolean isCNV(double value) {
        return value > 200d
    }
    
    int filteredTooSmall = 0
    
    int rawSignificantRegions = 0
    
    double minimumSignificantWidthFraction = 0.005
    
    /**
     * This method is called back by the thresholder each time a new contiguous region of
     * significant bases is observed (notionally, a candidate CNV).
     * <p>
     * For each such region, we expand the region to the cooresponding window size, centered around
     * the maximum of the signal. The expanded region is then combined with any other such regions already
     * detected that overlap it.
     * 
     * @param range     the range that was detected
     * @param width     the width of the corresponding detection window
     */
    void onNewCNVRegion(ThresholdRange range, int width) {
        RangeIndex cnvs = cnvs.get(width, new RangeIndex())
        
        ++rawSignificantRegions
        
        if(range.size()<minimumSignificantWidthFraction * width) {
            ++filteredTooSmall
            return
        }
        
        // Expand the call into a <width> wide region around the highest point
        GRange cnvPoint = new GRange(range.maxIndex, range.maxIndex, range)
        GRange cnvRange = cnvPoint.widen(width>>1i) // expand by width/2 on each side
        
        // Find all the overlapping calls: we will keep only the best
        Map<Boolean,List<GRange>> betterCnvs = (Map<Boolean,List<GRange>>)cnvs.getOverlaps(cnvRange.from, cnvRange.to).groupBy { IntRange r ->
            ThresholdRange otherRange = (ThresholdRange)((GRange)r).extra
            return otherRange.stats.max > range.stats.max
        }
        
        log.info "CNV ${range.from}-${range.to} (${range.size()}/${width}bp) has overlaps: better(${betterCnvs[true]?.size()}), worse(${betterCnvs[false]?.size()})"
        
        if(betterCnvs[true] == null)
            cnvs.add(cnvRange)
    }
    
    void updateThresholders(Contrast contrast, List<Double> deletionLRs) {
        if(thresholders.is(null)) {
            thresholders = []
            
            int index = 0
            for(Double lr : deletionLRs) {
                ++index
                
                int width = widths[index]
                Thresholder thresholder = 
                    new Thresholder(index:contrast.position)
                        .threshold(this)
                        .andThen((Closure<ThresholdRange>){ ThresholdRange tr ->
                            onNewCNVRegion(tr,width)
                        })
                thresholders.add(thresholder)
            }
        }
        
        for(int i=0; i<deletionLRs.size(); ++i) {
            thresholders[i].update(contrast.position,deletionLRs[i])
        }
    }

    @Override
    public Object include(double value, Object additional) {
        return isCNV(value);
    }
}
