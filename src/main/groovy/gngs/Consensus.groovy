package gngs
import groovy.transform.CompileStatic;

/**
 * Stores a set of base counts for a base read at a given position,
 * and supports scoring metrics.
 * 
 * @author simon
 */
@CompileStatic
class BaseCounts {
    
    BaseCounts(List<String> seqs, int pos) {
      for(int j=0;j<seqs.size();++j) {
        String seq = seqs[j]
        if(seq.size()<=pos)
            continue
        switch(seqs[j][pos]) {
          case 'A': ++a; break;
          case 'C': ++c; break;
          case 'T': ++t; break;
          case 'G': ++g; break;
          default: ++n;
        }
      }
    }
    
    int a=0,t=0,c=0,g=0,n=0
    
    String toString() { "cons: ${getConsensus()}, a:$a, t:$t, c:$c, g:$g, n:$n, score: ${getScore()}"  }
    
    char consensus = 'n'
    
    double getScore() {
      if(score < 0)
        score()
      return score
    }
    
    char getConsensus() {
      if(score < 0)
        score()
      return consensus
    }
     
    double score = -1.0d
    
    double score() {
      
      double total = Math.max(a+c+t+g,2) -1
      
      if(a>=t && a>=c && a>=g) {
        consensus='A'
        score = (a-1) / total
      }
      else
      if(c>=t && c>=g) {
        score = (c-1) / total
        consensus='C'
      }
      else
      if(g>=t) {
        score = (g-1) / total
        consensus='G'
      }
      else {
        consensus='T'
        score = (t-1) / total
      }
    }
    
    double bestScore() {
      Math.max(Math.max(Math.max(a,t),c),g)/(a+t+c+g)
    }
}

/**
 * Support to determine a consensus sequence from a set of sequences
 * 
 * @author simon
 */
@CompileStatic
class Consensus {
    
    BaseCounts [] baseCounts
    
    List<String> seqs
    
    int scanMax
    
    Consensus(List<String> sequences, int scanMax=Integer.MAX_VALUE) {
        
        this.seqs = sequences
        
        this.scanMax = Math.min(scanMax, seqs*.size().max())
        
        baseCounts = new BaseCounts[sequences.max { it.size() }.size()]
    }
    
    String bases = null
    
    double score = 0.0d
    
    Consensus build() {
        StringBuilder cons = new StringBuilder()
        for(int i=0; i<scanMax; ++i) {
          BaseCounts bc = new BaseCounts(seqs, i)
//          println "Base counts at $i: $bc"
          
          cons.append(bc.getConsensus())
          score += bc.score
        }
        bases = cons.toString()
        return this
    }

}
