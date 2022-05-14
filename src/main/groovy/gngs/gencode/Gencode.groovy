package gngs.gencode

import gngs.*
import graxxia.FileLike
import groovy.transform.CompileStatic

import htsjdk.tribble.index.Index
import htsjdk.tribble.readers.*

class Gencode {
    
    private final static int EQUALS = '=' as char

    Object source
    
    ProgressCounter progress = null
    
    Regions geneRegions = new Regions()
    
    Map<String, Gene> genes = new HashMap<String, Gene>(30000)
    
    Map<String,Feature> features = new HashMap<String, Gene>(30000)
    
    boolean loaded = false
    
    Index index
    
    List<Region> unresolvedEntriesAtPosition = []
    
    Gencode(def fileLike) {
        this.source = fileLike
    }
    
    @CompileStatic
    Gencode load() {
        
        try {
            List columnNames = [ 'id','source','type','start','end','score','strand','phase','attributes']

            new RangedData(FileLike.readerFactory(source), 0, 3, 4).load(columnNames: columnNames, commentChar: '#') { Region region ->
            
                if(!this.progress.is(null))
                    this.progress.count()
            
                parseEntry(region)
                
                return false
            } 
        }
        finally {
            if(!this.progress.is(null))
                this.progress.end()
            this.loaded = true
        }
        
        return this
    }
    
    @CompileStatic
    private Region parseEntry(final Region region) {
        
        if(!unresolvedEntriesAtPosition.isEmpty())  {
            if(region.from != unresolvedEntriesAtPosition[0].from) {
                
                List toResolve = unresolvedEntriesAtPosition
                unresolvedEntriesAtPosition = []
                
                // Resolve them
                for(Region unresolvedRegion in toResolve) {
                    addEntry(unresolvedRegion)
                }
                unresolvedEntriesAtPosition.clear()
            }
        }

        addEntry(region)
    }

    @CompileStatic
    private Region addEntry(final Region region) {
        Map<String,Object> attributes = parseAttributes(region)
        
        println "Parse $region: " + attributes
        if(region.from == 1335278) {
            println "Found it: $region"
        }
        
        final String type = region['type']
        Feature feature = null
        final String regionId = attributes['ID']
        if(type == 'gene') {
            region.properties.remove('attributes')
            this.genes[(String)attributes.gene_name] = feature = new Gene(region, regionId, (String)attributes.gene_name)
            this.geneRegions.addRegion(region)
        }
        else
        // coding exons appear both with the exon and CDS designation
        // UNLESS they are part of the UTR
        if(type == 'CDS' || type == 'three_prime_UTR' || type == 'five_prime_UTR') {
            feature = new Exon(region, regionId)
        }
        else
        if(type == 'transcript') {
            feature = new Transcript(region, regionId)
        }
        else {
            // This must be a non-exon / transcript / gene 
            return 
        }

        this.features[regionId] = feature

        addToParent(feature, attributes, region)

        return region
    }
    
    private final void addToParent(final Feature feature, final Map<String, Object> attributes, final Region region) {
        String parentId = attributes['Parent']
        if(parentId.is(null)) 
            return

        Feature parent = features[parentId]
        // if we are parsing a whole GFF3 file then parent should not be null, but in cases where we are parsing
        // just a subset of it, it can be
        if(!parent.is(null)) {
            feature.parent = parent
            parent.addChild(feature)
        }
        else {
            unresolvedEntriesAtPosition.add(region)
        }
    }

    private final static LinkedHashMap parseAttributes(Region region) {
        Map<String,Object> attributes = [:]
        String[] rawAttributes = ((String)region['attributes']).split(';')
        for(String rawAttribute in rawAttributes) {
            int equalsIndex = rawAttribute.indexOf(EQUALS)
            if(equalsIndex>=0) {
                attributes.put(rawAttribute.substring(0,equalsIndex), rawAttribute.substring(equalsIndex+1))
            }
            else {
                attributes.put(rawAttribute, Boolean.TRUE)
            }
        }
        return attributes
    }
    
    @CompileStatic
    void loadRegion(IRegion query) {
        
        if(!new File(this.source.toString() + '.tbi').canRead())
            throw new IllegalArgumentException("$source does not appear to be a readable tabix indexed file")
        
        TabixReader tbr = new TabixReader(this.source.toString(), this.source.toString() + '.tbi')
        TabixReader.Iterator i = tbr.query(query.chr, query.range.from, query.range.to)
        while(true) {
            List<String> line = i.next()?.tokenize('\t')
            if(line.is(null))
                break
                
            // NOTE: we use zero-based coordinates (as does BED format) but GFF3 uses 1 based
            // Therefore we have to subtract 1
            final int startPos = line[3].toInteger() -1
            
            // Even though end position is 1 based, because our range is
            // inclusive we do not subtract 1
            final int endPos = line[4].toInteger()

            Region region = new Region(line[0], startPos, endPos, type: line[2], source: line[1], attributes: line[8], strand: line[6])
            parseEntry(region)
        }
    }
    
    Gene getGeneRegion(String geneSymbol) {
        return genes[geneSymbol]
    }
}
