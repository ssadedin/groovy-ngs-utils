package gngs.gencode

import gngs.*
import graxxia.FileLike
import groovy.transform.CompileStatic

import htsjdk.tribble.index.Index
import htsjdk.tribble.readers.*
/**
 * Class implementing parsing of Gencode GFF3 format
 * <p>
 * This class passes a GFF3 file and creates a 3 level structure made up of
 * {@link Gene}, {@link @Transcript} and {@link Exon} objects, representing
 * the core elements within the GFF3. Sub elements of these within the file are
 * ignored.
 * <p>
 * There are two basic ways of using the class, and which one is better depends on the
 * kind of access you expect to make.
 * <ul>
 * <li>Fully parsing the file (requires significant memory, and takes quite a few seconds)
 * <li>Indexed access to look up specific genes on the fly
 * </u>
 * Full parsing is executed by simply creating the Gencode object and calling the {@link load()}
 * method to load it:
 * <pre>
 * def gencode = new Gencode()
 * 
 * </pre>
 * 
 * <p>To use indexed access, the gencode source file must be first sorted by position and then indexed using
 * tabix. A suitable command to do that is as follows:
 * <pre>
 * (zgrep ^"#" gencode.v40.basic.annotation.gff3.gz; zgrep -v ^"#" gencode.v40.basic.annotation.gff3.gz | sort -k1,1 -k4,4n) | bgzip > gencode.v40.basic.annotation.gff3.bgz
 * tabix -p gff gencode.v40.basic.annotation.gff3.bgz
 * </pre>
 * 
 * <p>After sorting and indexing, you may load specific regions by using the {@link #loadRegion} method:</p>
 * <pre>
 *     Gencode gencode = new Gencode(gencode)
 *     gencode.loadRegion(new Region('chr1:1335276-1349418'))
 * </pre>
 * 
 * @author simon.sadedin
 */
class Gencode implements GeneAnnotationSource {
    
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
       
        final String type = region['type']
        Feature feature = null
        final String regionId = attributes['ID']
        if(type == 'gene') {
            region.properties.remove('attributes')
            Gene gene = new Gene(region, regionId, (String)attributes.gene_name)
            if(attributes.hgnc_id)
                gene.hgnc_id = ((String)attributes.hgnc_id).split(':')[-1]
            gene.type = (String)attributes.gene_type
            
            feature = gene
            this.genes[(String)attributes.gene_name] = feature
            this.geneRegions.addRegion(region)
        }
        else
        // coding exons appear both with the exon and CDS designation
        // UNLESS they are part of the UTR
        if(type == 'CDS' || type == 'three_prime_UTR' || type == 'five_prime_UTR') {
            feature = new Exon(region, regionId)
            if(type == 'CDS')
                feature.coding = true
        }
        else
        if(type == 'transcript') {
            feature = new Transcript(region, regionId)
        }
        else {
            // This must be a non-exon / transcript / gene 
            return 
        }

        region['feature'] = feature

        this.features[regionId] = feature

        addToParent(feature, attributes, region)

        return region
    }
    
    @CompileStatic
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

    @CompileStatic
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
    
    @CompileStatic
    Gene getGeneFeature(String geneSymbol) {
        return genes[geneSymbol]
    }
    
    List<Gene> getGeneFeatures(IRegion region) {
         this.geneRegions.getOverlapRegions(region)*.feature       
    }
    
    @CompileStatic
    List<String> getGenes(IRegion region)  {
        List<Gene> genes = getGeneFeatures(region)
        return genes*.symbol
    }
    
    @CompileStatic
    Regions getExons(String geneSymbol, boolean codingOnly=true) {
        Gene gene = getGeneFeature(geneSymbol)
        List<IRegion> allRegions = (List<IRegion>)gene.children.collect { Feature transcript ->
            transcript*.children
        }.flatten()
        
        if(codingOnly)
            allRegions.removeIf { Feature exonFeature -> !((Exon)exonFeature).coding }

        Regions allExons = new Regions(allRegions)
        return allExons.reduce()
    }
    
    Map<String,Integer> getCDS(final IRegion region)  {
        List<Gene> genes = getGeneFeatures(region)
        return genes.collectEntries {gene -> 
            List<Exon> allCodingExons = gene.children*.children.flatten().grep { it.coding }
            Regions exonRegions = new Regions(allCodingExons)
            [
                gene.symbol,
                exonRegions.reduce().intersect(region)*.size().sum()?:0
            ]
        }
    }
}
