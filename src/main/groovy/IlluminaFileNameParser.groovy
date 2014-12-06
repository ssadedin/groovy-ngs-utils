import java.text.ParseException;


class IlluminaFileNameParser {
    
    /**
     * Mainly to allow some simple tests to work, we can let parsing be less strict
     * In that case, the sample name may be determined purely as all characters up to the 
     * first underscore. This is not satisfactory in a more strict environment.
     */
    boolean strict = Boolean.parseBoolean(System.getProperty("gngs.strict", "true"))
    
    static enum DIALECT { DEFAULT, MGHA }
    
    DIALECT dialect = DIALECT.DEFAULT
    
    /**
     * Value assigned to fields that could not be determined
     */
    public static final String UNKNOWN = "Unknown"
    

    /**
     * Examples of the various formats we try to parse:
     * 
     *  L10248_S12_L001_R2_001.fastq.gz
     *  GHS025_D1EH2ACXX_ATCACG_L001_R1.fastq.bam
     *  008_11_H14NGADXX_L1.1.fastq.trim.atrim.reorder.bam
     *
     * @param fileName
     * @return
     */
    Map parse(String fileName) { 
        
        def result = [
            sample: UNKNOWN,
            lane: UNKNOWN,
            unit: UNKNOWN
        ]
        
        // Don't include any relative portion of the path
        fileName = new File(fileName).name
        
        // Sometimes the index bar code is in the file name
        // We don't care about it, so it simplifies things to remove it
        def barcodeMatch = (fileName =~ /_([ATCG]{6,8})_/)
        if(barcodeMatch) {
            fileName = fileName.replaceAll("_"+barcodeMatch[0][1], "")
        }
        
        // If there is a run identifier, remove it
        fileName = fileName.replaceAll("_RUN[0-9*]_","_")
        
        def laneMatch = (fileName =~ /_(L[0-9]*)[._]/)
        if(laneMatch) {
            result.lane = laneMatch[0][1]
        }
        
        // Illumina machine names seem to end with XX, so make use of this to locate
        // the machine name part of the file name
        def machineNameMatch = (fileName =~ /_([A-Z0-9]*XX)_/)
        if(machineNameMatch) {
            // Sample name is everything up to the machine name
            String machineName = machineNameMatch[0][1]
            result.unit = result.lane + "." + machineName
            result.sample = fileName.substring(0,fileName.indexOf(machineName)-1)
        }
        else {
            if(result.lane != UNKNOWN) {
                
              // Dedicated parsing for MGHA anonymised samples which lose info from file names
              if((dialect == DIALECT.MGHA) && (fileName ==~ /S[0-9]{10}_[0-9]{2}_[a-z0-9]{10}.*/)) {
                  result.sample = fileName[0..10]
              }
              else 
                  result.sample = fileName.substring(0,fileName.indexOf(result.lane))
            }
            else {
                if(strict)
                    throw new RuntimeException("Unable to identify machine name or lane in file name $fileName")
                    
              result.sample = fileName.substring(0, fileName.indexOf("_"))
            }
        }
        
        result.sample = result.sample.replaceAll('_$', '')
        
        return result
    }
}
