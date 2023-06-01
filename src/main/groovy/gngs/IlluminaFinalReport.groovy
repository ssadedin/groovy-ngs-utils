package gngs

import groovy.transform.CompileStatic

@CompileStatic
class IlluminaSNPData implements IRegion {
    String snpName
    String chr
    int pos

    double gtScore
    double gcScore

    /**
     * Call of first allele relative to forward strand
     */
    String allele1


    /**
     * Call of first second allele relative to forward strand
     */
     String allele2

     double logR

     double theta

     double R

    @Override
    public IntRange getRange() {
        return pos..(pos+60) // hack - should use probe direction and real length of probe 
    }
}


class IlluminaFinalReport {
    
    List<IlluminaSNPData> data = new ArrayList(500000)
    
    String toString() {
        "IlluminaFinalReport(${data.size()} probes)"
    }

    static IlluminaFinalReport parse(path) {
        Utils.reader(path) {
            parseReader(it)
        }
    }

    @CompileStatic
    static IlluminaFinalReport parseReader(Reader r) {
        
        String line = null
        try {
            // File constists of two sections:
            // [Header]
            // ... headers
            // [Data]
            // ... data rows
            
            IlluminaFinalReport ifr = new IlluminaFinalReport()

            String headerLine = r.readLine()
            if(headerLine != "[Header]")
                throw new IllegalArgumentException('File expected to start with [Data] element')

            // First read all the headers by reading lines until we get to a [Data]
            List<String> headers = []
            while((line = r.readLine()) != "[Data]") {
                if(line.is(null))
                    break
                headers << line
            }

            if(line != "[Data]")
                throw new IllegalArgumentException('Unable to locate [Data] section in file')

            String columnHeaderLine = r.readLine()
            final List<String> hds = columnHeaderLine.split("\t") as List
            
            // Chr  Position    SNP Name    GT Score    GC Score    Allele1 - AB    Allele2 - AB    Allele1 - Top   Allele2 - Top   Allele1 - Forward   
            // Allele2 - Forward   Allele1 - Design    Allele2 - Design    Theta   R   X Raw   Y Raw   X   Y   Log R Ratio B Allele Freq   SNP Aux SNP ILMN Strand Top Genomic Sequence    Customer Strand

            int snpNamePos = hds.indexOf('SNP Name')
            int chrPos = hds.indexOf('Chr')
            int posPos = hds.indexOf('Position')
            int gtScorePos = hds.indexOf('GT Score')
            int gcScorePos = hds.indexOf('GC Score')
            int allele1Pos = hds.indexOf('Allele1 - Forward')
            int allele2Pos = hds.indexOf('Allele2 - Forward')
            int thetaPos = hds.indexOf('Theta')
            int rPos = hds.indexOf('R')
            int logRPos = hds.indexOf('Log R Ratio')

            // Finally, read the data
            while((line = r.readLine()) != null) {
                String [] fields = line.split('\t')
                IlluminaSNPData data = new IlluminaSNPData()
                data.chr = 'chr' + fields[chrPos]
                data.pos = fields[posPos].toInteger()
                data.allele1 = fields[allele1Pos]
                data.allele2 = fields[allele2Pos]
                data.gtScore = fields[gtScorePos].toDouble()
                data.gcScore = fields[gcScorePos].toDouble()
                data.theta = fields[thetaPos].toDouble()
                data.logR = fields[logRPos].isEmpty() ? Double.NaN : fields[logRPos].toDouble()
                data.R = fields[rPos].toDouble()
                ifr.data.add(data)
            }
            
            return ifr
        }
        catch(Exception e) {
            if(line != null) {
                throw new IOException("Failed to read IlluminaFinalReport parsing line:\n\n$line", e)
            }
            else
                throw e
        }
    }
}
