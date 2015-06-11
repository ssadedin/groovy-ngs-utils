import graxxia.*

class BEDToolsCoverage {
    
    Regions covs = new Regions()
    
    IntegerStats allStats = new IntegerStats(1000)
    
    IntegerStats inStats = new IntegerStats(1000)
    
    BEDToolsCoverage(String file) {
        this(file,null)
    }
    
    BEDToolsCoverage(String file, Regions mask) {
        
        // Sniff the file to see if there are 5 columns or 6 - if the original bed file had an id 
        // column then the BedTools output also has it
        def firstLine = new File(file).newReader().readLine().split('\t')
        
        List<String> cols = ["chr","start","end","offset","cov"]
        boolean hasId = false
        if(firstLine.size()>5) {
            hasId = true
            cols = ["chr","start","end","id","offset","cov"]
        }
        
        TSV covFile = new TSV(file, columnNames:cols, readFirstLine:true)
        Region r = null
        List regionCovs = []
        boolean inMask = false
        for(cov in covFile) {
           if(r == null || r.chr != cov.chr || r.from != cov.start || r.to != cov.end) { // new region
               if((r != null) && inMask) // add previous region
                   covs.addRegion(r)
               r = new Region(cov.chr, cov.start..cov.end)
               r.covs = []
               
               inMask = (mask == null || r.overlaps(mask))
               
               if(hasId)
                   r.id = cov.id
           }
           if(inMask)
               inStats.addValue(cov.cov)
           allStats.addValue(cov.cov)    
           r.covs.add(cov.cov)
           
        }
        if((r != null) && inMask) // add previous region
            covs.addRegion(r)
    }
    
    List<Integer> coverage(Region query) {
        List<IntRange> overlaps = covs.getOverlaps(query)
        List<Integer> values = overlaps.sum { IntRange range ->
          int startOffset = query.from > range.from ? query.from - range.from : 0
          int endOffset = query.to < range.to ? range.to - query.to : 0
          
          if(startOffset == 0 && endOffset == 0)
              range.extra.covs
          else
              range.extra.covs.subList(startOffset, range.extra.covs.size() - endOffset)
        }
        
        return (values == null ? [] : values)
    }
    
    CoverageStats stats(Region query) {
        return new CoverageStats(5000,coverage(query))
    }
}
