package gngs

import groovy.transform.CompileStatic

/**
 * Static constant information about various genomes used to enable auto-sniffing etc
 * 
 * @author simon.sadedin
 */
@CompileStatic
class GenomeInfo {
        static Map hgMap = [
            247249719 : "hg18",
            249250621 : "hg19",
            248956422 : "hg38"
        ]
        
        static Map grcMap = [
            249250621 : "GRCh37",
            248956422 : "GRCh38"
        ]
        
        static Map mouseMap = [
            195471971 : "mm10",
            197195432 : "mm9",
            197069962 : "mm8",
            194923535 : "mm7"
        ]
        
        /**
         * @param contigSizes
         * @return 38, 37 or 18 indicating the coordinate system of the genome version 
         *         for the given contig sizes
         */
        static Integer humanGenomeCoordinateVersion(final Map<String, Integer> contigSizes) {
            String build = sniffGenomeBuild(contigSizes)
            if(build in ["hg38", "GRCh38"])
                return 38

            if(build in ["hg19", "GRCh37"])
                return 37

            if(build in ["hg18"])
                return 36

            return -1
        }

        /**
         * Infers the identity of the genome build used based on given contig sizes.
         * 
         * @param contigSizes
         * @return  String identifier for genome build
         */
        static String sniffGenomeBuild(final Map<String, Integer> contigSizes) {
            contigSizes.any { it.key.startsWith('chr') } ? 
                hgMap.getOrDefault(contigSizes['chr1'], mouseMap[contigSizes['chr1']])
            :
                grcMap[contigSizes['1']]
        }
}
