/*
 *  Groovy NGS Utils - Some simple utilites for processing Next Generation Sequencing data.
 *
 *  Copyright (C) 2014 Simon Sadedin, ssadedin<at>gmail.com
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
public class Pileup implements Iterator<Pileup> {
    
    public String chr;
    
    public int start;
    
    public int end;
    
    public int position;
    
    public List<PileupState> alignments = new ArrayList<PileupState>(1000); 
    
    private SAMRecordIterator readIterator;
    
    private SAMRecord nextRead;
    
    private SAMFileReader samFile;
    
    public Pileup() {}
    
    public Pileup(SAMFileReader samFile, String chr, int start, int end) {
        this.start = start;
        this.end = end;
        this.chr = chr;
        this.readIterator = samFile.query(chr, start, end, false);
        this.samFile = samFile;
        this.position = -1;
    } 
    
    public Pileup next() {
        
        if(!hasNext()) {
            this.alignments = new ArrayList<PileupState>();
            return this;
        }
        
         do {
          if(position != -1)
              ++position;
          
          SAMRecord r = nextRead;
          while(true) {

              // Skip reads until we find a valid one
              while(r == null && this.readIterator.hasNext()) {
                  r = (SAMRecord)readIterator.next();
                  
                  // Ignore reads that are unmapped or not mapped uniquely   
                  if(r != null && (r.getReadUnmappedFlag() || r.getMappingQuality()==0))
                      r = null;
              }

              if(position == -1) {
                  if(r == null) // no valid reads over the current position
                      break;
                  position = r.getAlignmentStart();
              }

              // Pile up all the reads at this position
              if(r != null && r.getAlignmentStart() == position) {
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
              // System.out.println(this.position + " : " + (char)s.base);
          }
          this.alignments.removeAll(remove);
          this.nextRead = r;        
        } 
        while(position > 0 && position < start) ;
        
        return this;
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
    
    /**
     * Return the count of reads supporting the specified base at the current
     * position.
     * 
     * @param base	ascii int code for A,T,C,G,N, etc.
     * @return		count of reads supporting base
     */
    public int countOf(byte base) {
        int count = 0;
        for(PileupState s : alignments) {
            if(s.base == base)
                ++count;
        }
        return count;
    }

    @Override
    public boolean hasNext() {
        return position < end && !(!readIterator.hasNext() && position == -1);
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException();
    }
    
    public void close() {
        this.readIterator.close();
        this.samFile.close();
    }
}
