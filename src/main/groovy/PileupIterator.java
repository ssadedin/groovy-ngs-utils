import groovy.lang.Closure;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

/**
 * Overall class containing pileup state
 * 
 * @author simon.sadedin@mcri.edu.au
 */
public class PileupIterator implements Iterator<PileupIterator.Pileup> {
    
    /**
     * Inner class that captures pileup state at a specific position
     * 
     * @author simon.sadedin@mcri.edu.au
     */
    public class Pileup {
        
        public int position = -1;
        
        Pileup(int position) {
            this.position = position;
        }
        
        /**
         * Returns a summary of the coverage at the location:
         * 
         * [ A, T , C, G, Deletion, total ]
         * 
         * @return
         */
        public int[] getSummary() {
            int a=0,t=0,c=0,g=0,d=0;
            
            for (PileupState s : alignments) {
                switch(s.base) {
                
                  case 'A':
                      ++a;
                      break;
                  case 'T':
                      ++t;
                      break;
                  case 'C':
                      ++c;
                      break;
                  case 'G':
                      ++g;
                      break;
                  case 0:
                      ++d;
                      break;
                }
            }
            return new int[] { a,t,c,g,d,a+t+c+g+d};
        }
        
        /**
         * Convenience method to return counts of each base as a map indexed
         * by the base as a string. This is not efficient, but useful in non
         * performance sensitive contexts.
         * @return  map of count indexed by base
         */
        public Map<String, Integer> getSummaryAsMap() {
            int [] summary = getSummary();
            HashMap<String, Integer> result = new HashMap<String,Integer>();
            result.put("A", summary[0]);
            result.put("T", summary[1]);
            result.put("C", summary[2]);
            result.put("G", summary[3]);
            result.put("D", summary[4]);
            result.put("total", summary[5]);
            return result;
        }
        
        public int countOf(String baseString) {
            return countOf((byte)baseString.charAt(0));
        }
        
        public int countOf(byte base) {
            int count = 0;
            for(PileupState s : alignments) {
                if(s.base == base)
                    ++count;
            }
            return count;
        }
        
        public int countOf(Closure<Object> c) {
            int count = 0;
            for(PileupState s : alignments) {
                Object result = c.call(s);
                if(result != null && result == Boolean.TRUE)
                    ++count;
            }
            return count;
        } 
        
        public List<PileupState> getAlignments() {
            return alignments;
        }
    }
    
    public String chr;
    
    public int start;
    
    public int end;
    
    public int iteratorPosition;
    
    public List<PileupState> alignments = new ArrayList<PileupState>(1000); 
    
    private SAMRecordIterator readIterator;
    
    private SAMRecord nextRead;
    
    private SAMFileReader samFile;
    
    private int minMappingQuality = 1;
    

    public PileupIterator() {}
    
    public PileupIterator(SAMFileReader samFile, String chr, int start, int end) {
        this.start = start;
        this.end = end;
        this.chr = chr;
        this.readIterator = samFile.query(chr, start, end, false);
        this.samFile = samFile;
        this.iteratorPosition = -1;
    } 
    
    public Pileup next() {
        
        if(!hasNext()) {
            this.alignments = new ArrayList<PileupState>();
            return new Pileup(this.iteratorPosition);
        }
        
         do {
          if(iteratorPosition != -1)
              ++iteratorPosition;
          
          SAMRecord r = nextRead;
          while(true) {

              // Skip reads until we find a valid one
              while(r == null && this.readIterator.hasNext()) {
                  r = (SAMRecord)readIterator.next();
                  
                  // Ignore reads that are unmapped or mapped below the threshold set
                  if(r != null && (r.getReadUnmappedFlag() || r.getMappingQuality()<this.minMappingQuality))
                      r = null;
              }

              if(iteratorPosition == -1) {
                  if(r == null) // no valid reads over the current position
                      break;
                  iteratorPosition = r.getAlignmentStart();
              }

              // Pile up all the reads at this position
              if(r != null && r.getAlignmentStart() == iteratorPosition) {
                  alignments.add(new PileupState(r));
                  nextRead = null;
                  r = null;
                  continue;
              }
              else 
                  break;
          }

          // The next read is past the current position - we should return
          // First though, remove all the reads that do not overlap this position
          List<PileupState> remove = new ArrayList<PileupState>(this.alignments.size());
          for(PileupState s :  this.alignments) {
              if(!s.next()) // next returns false if the position doesn't overlap
                  remove.add(s);
          }
          this.alignments.removeAll(remove);
          this.nextRead = r;        
        } 
        while(iteratorPosition > 0 && iteratorPosition < start) ;
        
        return new Pileup(this.iteratorPosition);
     }
    

    @Override
    public boolean hasNext() {
        return iteratorPosition < end && !(!readIterator.hasNext() && iteratorPosition == -1);
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException();
    }
    
    public void close() {
        this.readIterator.close();
        this.samFile.close();
    }
    
    public int getMinMappingQuality() {
        return minMappingQuality;
    }

    public void setMinMappingQuality(int minMappingQuality) {
        this.minMappingQuality = minMappingQuality;
    }
}
