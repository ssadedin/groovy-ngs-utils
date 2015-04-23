
class BEDToolsCoverage {
    
    Regions covs = new Regions()
    
    BEDToolsCoverage(String file) {
        TSV covFile = new TSV(file, columnNames:["chr","start","end","gene","offset","cov"], readFirstLine:true)
        Region r = null
        List regionCovs = []
        for(cov in covFile) {
           if(r == null || r.chr != cov.chr || r.from != cov.start || r.to != cov.end) { // new region
               if(r != null) // add previous region
                   covs.addRegion(r)
               r = new Region(cov.chr, cov.start..cov.end)
               r.covs = []
               r.gene = cov.gene
           }
           r.covs.add(cov.cov)
        }
        if(r != null) // add previous region
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
