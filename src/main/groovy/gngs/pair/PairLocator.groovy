/*
 *  Groovy NGS Utils - Some simple utilites for processing Next Generation Sequencing data.
 *
 *  Copyright (C) 2018 Simon Sadedin, ssadedin<at>gmail.com
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
package gngs.pair

import java.io.Writer
import java.util.concurrent.ConcurrentHashMap
import gngs.*

import groovy.transform.CompileStatic
import groovy.util.logging.Log
import groovyx.gpars.actor.Actor
import groovyx.gpars.actor.DefaultActor
import htsjdk.samtools.SAMRecord

@CompileStatic
class Paired<PairType> {
	
	final static byte R1_READ_NEGATIVE_STRAND_FLAG = 0x1
	final static byte R2_READ_NEGATIVE_STRAND_FLAG = 0x2
	final static byte REORDER_IN_OUTPUT = 0x4
	
	public String key
    public PairType r1
    public PairType r2
	int flags = 0

    public Paired(String key, PairType r1, PairType r2, boolean r1NegStrand, boolean r2NegStrand, boolean flipOrder, int rand) {
        super();
		this.key = key
        this.r1 = r1;
        this.r2 = r2;
		this.flags |= (r1NegStrand ? R1_READ_NEGATIVE_STRAND_FLAG : 0)
		this.flags |= (r2NegStrand ? R2_READ_NEGATIVE_STRAND_FLAG : 0)
		this.flags |= (flipOrder ? REORDER_IN_OUTPUT : 0)
        this.flags |= (rand & 0xFFF8)
    }
	
	boolean getR1NegativeStrandFlag() {
		flags & R1_READ_NEGATIVE_STRAND_FLAG
	}
	
	boolean getR2NegativeStrandFlag() {
		flags & R2_READ_NEGATIVE_STRAND_FLAG
	}
    
    boolean getReorderRequired() {
        flags & REORDER_IN_OUTPUT
    }
}

@CompileStatic
@Log
class PairLocator<PairType extends ReadPair> extends RegulatingActor<List<SAMRecord>> {
    
    /**
     * Read pairs we are indexing, keyed on name
     */
    Map<String,PairType> buffer 
    
    RegulatingActor consumer
    
    int received = 0
    
    int paired
    
    int chimeric
    
    Regions regions
    
    String debugRead = null // "SN7001291:342:HFMC7BCXX:2:1113:1980:38527"
    
    Set<Integer> chromosomesWithReads
    
    boolean compact = true
    
    /**
     * The tag from which to extract base quality scores (use actual base qualities if null)
     */
    String baseQualityTag = null

    /**
     * Random number generator used to randomize read output order
     */
    Random rand = new Random(0)
        
    PairLocator(RegulatingActor<Paired> consumer, Set<Integer> chromosomesWithReads) {
        super(20000,50000)
        this.consumer = consumer
        this.buffer = new HashMap(200000)
        this.chromosomesWithReads = chromosomesWithReads
    }
    
    void process(final List<SAMRecord> records) {
        for(SAMRecord r : records) {
            processRead(r)
        }
    }
    
    @Override
    public void onEnd() {
        process(assigned)
    }

    void processRead(final SAMRecord record) {
        
        ++received
        
        if(record.secondaryOrSupplementary)
            return
        
        final String readName = record.readName
        PairType pair = buffer.remove(readName)
        if(!pair.is(null)) {
            emitRead(pair, record)
            return
        }
        
        // If we will never mate this read, let's not store it
        if(!chromosomesWithReads.is(null) && !chromosomesWithReads.contains(record.mateReferenceIndex))
            return
        
        pair = new SAMRecordPair(r1:record)
        
        // Is it or its mate located in the desired regions?
        if(notInRegion(pair, readName == debugRead)) {
            if(readName == debugRead)
                log.info "$readName filtered out"
            return
        }
        
        if(debugRead != null && (readName == debugRead)) {
            log.info "Buffer: $record"
        }
            
        // Pair not found - create a new one
        // Note: we avoid creating the compact read pair for reads we will not actually buffer because
        //       the constructor may do expensive compression which is wasted if we then throw the read 
        //       away 
        if(compact)
            pair = new CompactReadPair(record, this.baseQualityTag)
       
        if(pair.chimeric)
            ++this.chimeric
            
        buffer.put(readName,pair)
    }

    private emitRead(final ReadPair pair, final SAMRecord record) {
        final String readName = record.readName
        if(!debugRead.is(null) && (readName == debugRead)) {
            log.info "Paired: $record"
        }
        if(pair instanceof SAMRecordPair)
            pair.r2 = record

		Paired p =new Paired(
    			record.readName, 
    			pair, 
    			new CompactReadPair(record, this.baseQualityTag), 
    			record.mateNegativeStrandFlag, 
    			record.readNegativeStrandFlag,
                record.firstOfPairFlag,
                rand.nextInt()
            )

        consumer.sendTo(p)
        buffer.remove(readName)
        paired += 2
    }
    
    boolean notInRegion(PairType pair, boolean debug) {
        
        if(pair.unmapped) {
            return false
        }
        
        if(regions.is(null))
            return false
        
        return pair.notInRegions(regions)
    }
    
    List<SAMRecord> assigned = []
    
    /**
     * Assigns the given read to this locator to process.
     * <p>
     * The read will be buffered until there are 10 reads to process, and then passed to the
     * actor thread.
     * 
     * @param read
     */
    void assign(SAMRecord read) {
        assigned.add(read)
        if(assigned.size()>20) {
            this.send(new AcknowledgeableMessage(assigned, this.pendingMessageCount))
            this.assigned = new ArrayList(20)
        }
    }
    
    String toString() {
        'Locator: ' + [received,paired,buffer.size()].collect { Utils.human(it) }.join(',')
    }
}