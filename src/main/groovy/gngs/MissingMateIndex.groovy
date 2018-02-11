package gngs

import groovy.transform.CompileStatic
import groovy.util.logging.Log
import htsjdk.samtools.SAMRecord


class ScannableRegion {
    
    public Region region
    
    SortedMap<Long,List<SAMRecordPair>> matePositions 
}

@Log
class RegionScan {
    
    ScannableRegion region 
    List<SAMRecordPair> ancilliaryPairs  
    
    void markUnpaired() {
        int unpairedCount = 0
        for(SAMRecordPair ancilliaryPair in ancilliaryPairs) {
            if(!ancilliaryPair.hasBothReads()) {
                ancilliaryPair.flags.unpaired = true
                ++unpairedCount
            }
        }
        if(unpairedCount>20) // just an arbitrary debug threshold: some can be unpaired, but a lot is suspicious
            log.info "${unpairedCount}/${ancilliaryPairs.size()} reads were not paired after expected resolution"

    }
}

@CompileStatic
@Log
class MissingMateIndex {
    
    /**
     * A list of unresolved read pairs (that is, where one pair is missing)
     * ordered by the priority of resolving them.
     * <p>
     * We assume that genome is being scanned in sorted order, and thus
     * the priority is to resolve the "left most" or lowest genomic
     * position first, as this is most likely to hold up scanning.
     */
    PriorityQueue<SAMRecordPair> resolveQueue = new PriorityQueue(100000)
    
    /**
     * Counts of read pairs we are indexing, keyed on missing mate position
     */
    SortedMap<Long,List<SAMRecordPair>> matePositions = Collections.synchronizedSortedMap(new TreeMap())
    
    /**
     * Read pairs we are indexing, keyed on name
     */
    Map<String,SAMRecordPair> readIndex 
    
    /**
     * Index of ALL chimeric reads
     * <p>
     * Note: this index is not cleared by the {@link#clear} method
     * because it needs to persist between chromosomes.
     */
    Map<String,SAMRecordPair> chimericReadIndex 
    
    boolean verbose = false
    
    MissingMateIndex(int initialSize) {
        this.readIndex = new HashMap(initialSize)
        this.chimericReadIndex = new HashMap(initialSize)
    }
    
    /**
     * Index this mate
     * 
     * @param pair
     */
    synchronized void add(SAMRecordPair pair) {
        
       // Chimeric reads are indexed separately because they need
       // to persist between chromosomes to be resolved
       if(pair.chimeric)
           chimericReadIndex[pair.readName] = pair
       else
           readIndex[pair.readName] = pair
       
       String chr = pair.nonMissingRead.mateReferenceName
       if(chr == "*") {
           return // don't really know what to do with this
       }
       
       if(Region.isMinorContig(chr))
           return
           
       long xPos = XPos.computePos(chr, pair.missingReadPosition)
       matePositions.get(xPos,[]).add(pair)
    }
    
    synchronized SAMRecordPair getPair(String readName) {
        readIndex[readName]?:chimericReadIndex[readName]
    }
    
    synchronized void resolve(int pos, SAMRecordPair pair) {
        String chr = pair.nonMissingRead.mateReferenceName
        if(chr!="*") {
            if(!Region.isMinorContig(chr)) {
                long xPos = XPos.computePos(chr, pos)
                if(xPos in matePositions) {
                    List pairsAtPos = matePositions[xPos]
                    if(pairsAtPos.size()==1)
                        matePositions.remove(xPos)
                    else
                        pairsAtPos.remove(pair)
                }
            }
        }
        if(pair.chimeric)
            chimericReadIndex.remove(pair.readName)
        else
            readIndex.remove(pair.readName)
    }
    
    /**
     * Search for a region that is guaranteed to contain the missing read from the specified pair,
     * while also containing at least <code>minimumReads</code> reads that also need resolution.
     * 
     * @param pair              the read pair to search for the missing read for
     * @param regionScanSize    how far to look upstream and downstream of the missing read
     * @param minimumReads      the minimum number of reads to consider worthy of a scannable region
     * 
     * @return  A region to scan or null if no region satisfying the requirements exists
     */
    synchronized ScannableRegion searchForScannableRegion(final SAMRecordPair pair, final int regionScanSize, final int minimumReads) {
        
        if(Region.isMinorContig(pair.missingReadReferenceName))
            return null
        
        int readPos = pair.missingReadPosition
        long mateXPos = XPos.computePos(pair.missingReadReferenceName, readPos)
        SAMRecord record = pair.nonMissingRead
        SortedMap<Long,List<SAMRecordPair>> positions  
        
        synchronized(this) {
            positions = matePositions.subMap(mateXPos-regionScanSize, mateXPos+regionScanSize)
            if(positions.isEmpty()) {
                if(verbose)
                    log.info "No reads identified between ${mateXPos-regionScanSize}-${mateXPos+regionScanSize} despite $pair with contained mate (mate index size = ${matePositions.size()})"
                return null
            }
            else {
                if(verbose)
                    log.info "Identified ${positions.size()} resolvable read positions between ${mateXPos-regionScanSize}-${mateXPos+regionScanSize} for $pair (mate index size = ${matePositions.size()})"
            }
        }
        
        
        Region startRegion = XPos.parsePos(positions.firstKey())
        
        int start = startRegion.from
        
        assert start != 0 : "Search of ${matePositions.size()} positions near ${pair.missingReadPosition} with ${positions.size()} needed records should not have returned zero start position"
        
        if(start > readPos)
            start = readPos-1
        
        int end = startRegion.from  + (int)(positions.lastKey() - mateXPos)
        if(end < readPos)
            end = readPos + 1
        
        int closeReadCount = (int)(positions*.value*.size().sum()?:0i)
        if(closeReadCount>minimumReads) {
            return new ScannableRegion(
                    region:new Region(record.mateReferenceName, startRegion.from,  end), 
                    matePositions:positions)
        }
        else {
            return null
        }
    }
    
    synchronized RegionScan allocateScannableRegion(SAMRecordPair pair, int regionScanSize, int minPositions) {
        
        RegionScan scan = new RegionScan()
        
        scan.region = this.searchForScannableRegion(pair, regionScanSize, minPositions) 
        if(scan.region == null) { 
             return null
         }
                 
         // Since we will find them, remove all the pairs from the queue
         // so we don't waste time processing them again later
         int oldSize = this.matePositions.size()
         int countRemoved = 0
         Iterator<Map.Entry<Long,List<SAMRecordPair>>> i = scan.region.matePositions.iterator()
                 
         // note, could be more than 1 read at a given position so size is >= to our estimate
         List pairs = new ArrayList((int)scan.region.matePositions.size()) 
         while(i.hasNext()) {
             List<SAMRecordPair> readPairs = i.next().value
             assert readPairs.size() > 0
             pairs.addAll(readPairs)
             i.remove()
             ++countRemoved
         }
         
         scan.ancilliaryPairs = pairs
         if(verbose)
             log.info "Removed $countRemoved mate positions ($oldSize => ${matePositions.size()})"
         return scan
    }
    
    synchronized int getResolveQueueSize() {
        this.resolveQueue.size()
    }
    
    synchronized void clear() {
        readIndex.clear()
        matePositions.clear()
    }
}
