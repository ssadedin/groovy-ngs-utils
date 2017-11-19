import gngs.Region
import gngs.Regions

/**
 * Reads output of GATK DepthOfCoverage interval summary and
 * makes it available as a Regions object
 * 
 * @author simon
 */
class GATKIntervalSummary extends Regions {
    
    File file
    
    static List<String> columns = [
        "target",
        "total_coverage",
        "average_coverage",
        "total_cvg",
        "mean_cvg",
        "granular_Q1",
        "granular_median",
        "granular_Q3",
        "perc_above_15"
    ]
    
    GATKIntervalSummary(String fileName) {
        this.file = new File(fileName)
    }
    
    GATKIntervalSummary load() {
        graxxia.TSV tsv = new graxxia.TSV(file.absolutePath, columnNames: columns)
        for(line in tsv) {
           if(!line.target.contains("-"))
               continue
               
           Region region = new Region(line.target) 
           for(prop in columns) {
               if(prop != "target")
                   region.setProperty(prop, line[prop].toDouble())
           }
           
           this.addRegion(region)
        }
        return this
    }
    
    Iterator<Region> iterator() {
      super.iterator()  
    }
}
