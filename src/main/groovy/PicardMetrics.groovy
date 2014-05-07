/**
 * A simple class for parsing information out of the metrics file written by Picard 
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class PicardMetrics {

    static Map<String,Object> parse(String fileName) {
        List<String> lines = new File(fileName).readLines()
        int index = lines.findIndexOf { it.startsWith("LIBRARY") }
        if(index < 0)
            throw new IllegalStateException("Unable to locate LIBRARY line in Picard metrics file ${fileName}")
            
        [ lines[index].split('\t')[1..-1], lines[index+1].split('\t')[1..-1]*.toFloat() ].transpose().collectEntries() 
    }
}
