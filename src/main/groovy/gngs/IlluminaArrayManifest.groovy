package gngs

import groovy.transform.CompileStatic
import groovy.transform.ToString

/**
 * An individual probe in an {@link gngs.IlluminaArrayManifest}.
 * 
 * @author simon.sadedin
 */
@CompileStatic
@ToString(excludes=['AlleleA_ProbeSeq'])
class IlluminaProbe implements IRegion {
    
    String IlmnID
    String Name
    String IlmnStrand
    String SNP
    String AlleleA_ProbeSeq
    String Chr
    int MapInfo
    String RefStrand
    boolean Intensity_Only
    byte type
    
    String getChr() { Chr }
    IntRange getRange() { RefStrand == '-' ?  MapInfo..<(MapInfo+AlleleA_ProbeSeq.size()) :  (MapInfo-AlleleA_ProbeSeq.size())..MapInfo}
}


/**
 * Support for parsing an Illumina Array Manifest file.
 * 
 * Each probe is parsed to a {@link IlluminaProbe} object that contains
 * most of the original fields from the file, but implements the 
 * {@link gngs.IRegion} interface, allowing it to support range operations.
 * <p> 
 * Probes may be accessed either via lookup from their name in the {@link #probes}
 * Map attribute, or they are indexed in a {@link gngs.Regions} object to search
 * by region (each {@link gngs.Region} object having a <code>probe</code> attribute
 * to allow full details of the probe to be retrieved.
 * 
 * @author simon.sadedin
 */
class IlluminaArrayManifest {
    
    Map<String,IlluminaProbe> probes
    
    Regions regions
    
    static IlluminaArrayManifest parse(Object fileLike) {
        Utils.reader(fileLike).withReader { r ->
            parseReader(r)
        }
    }

    @CompileStatic
    static IlluminaArrayManifest parseReader(Reader r) {
        
        // File starts with header like
        // Illumina, Inc.
        // [Heading]
        // Descriptor File Name,GDACyto_20047166_A2.bpm
        // Assay Format,Infinium LCG
        // Date Manufactured,4/30/2021
        // Loci Count ,2065511
        // [Assay]
        
        String line = ""
        int lociCount = -1
        int lineCount = 0
        while(line != null && !line.startsWith("[Assay]")) {
            line = r.readLine()
            if(line.startsWith("Loci Count")) {
                lociCount = line.tokenize(',')[-1].toInteger()
            }
            ++lineCount
        }
        
        if(line == null)
            throw new IllegalStateException("Manifest file ended prematurely : truncated file?")
            
        line = r.readLine()
        ++lineCount
        
        Map<String,Integer> headers = line.tokenize(',').indexed().collectEntries { [it.value.toLowerCase(), it.key] }
        if(!headers.containsKey('ilmnid'))
            throw new IllegalStateException("Expected header IlmnID not found in array manifest header line: $line. Please check file format is correct")
        
        Regions regions = new Regions()
        
        Map<String,IlluminaProbe> probes = lociCount > 0 ? new LinkedHashMap(lociCount+1) : [:]

        try {
            while((line = r.readLine()) != null) {
                
                if(line == "[Controls]")
                    break
                
                ++lineCount
                String [] fields = line.split(',')
                
                if(fields.size() < 15)
                    continue
                    
                final String chr = fields[headers.chr]
                IlluminaProbe probe = new IlluminaProbe(
                    IlmnID: fields[headers.ilmnid],
                    Name: fields[headers.name],
                    SNP: headers.snp != null ? fields[headers.snp] : null, 
                    AlleleA_ProbeSeq : fields[headers.allelea_probeseq],
                    Chr : chr.startsWith('chr') ? chr : 'chr'+chr,
                    MapInfo: fields[headers.mapinfo].toInteger(),
                    RefStrand: headers.refstrand ? fields[headers.refstrand] : (fields[headers.strand_fr] == "F" ? "+" : "-"),
                    Intensity_Only : headers.intensity_only ? fields[headers.intensity_only] != "0" : null,
                    type : fields[headers.addressb_id] ? 1 : 2
                )
                Region region = new Region(probe.chr, probe.range, probe:probe)
                regions.addRegion(region)
                
                probes[probe.Name] = probe
            }
        }
        catch(Exception e) {
            throw new ParseException(e, lineCount)
        }
        
        return new IlluminaArrayManifest(probes:probes, regions: regions)
    }
}
