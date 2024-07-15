package gngs.gencode

import com.google.common.collect.ImmutableList
import gngs.GRange
import gngs.IRegion
import groovy.transform.CompileStatic

@CompileStatic
class Feature<CHILD_TYPE extends Feature> implements IRegion {
    String chr
    IntRange range
    String id
    Feature parent
    List<CHILD_TYPE> children
    
    Feature(IRegion region, String id, Feature parent) {
        this.chr = region.chr
        this.range = new GRange(region.range.from, region.range.to, null)
        this.id = id
        this.parent = parent
    }
    
    void addChild(Feature child) {
        if(this.children.is(null)) {
           this.children = [] 
        }
        this.children.add(child)
    }
    
    List<Feature> getChildren() {
        if(!this.children || this.children.isEmpty())
            return Collections.emptyList()
        return this.children
    }
    
    int size() {
        return range.size()
    }
    
    String toString() {
        "${this.class.name.replace('gngs.gencode.','')} at $chr:$range.from-$range.to"
    }
}

@CompileStatic
class Gene extends Feature<Transcript> {
    String symbol
    String hgnc_id
    String type
    char strand
    
//    Transcript transcript
    public Gene(IRegion region, String id, String symbol, char strand) {
        super(region, id, null);
        this.symbol = symbol
        this.strand = strand
    }
}

@CompileStatic
class Transcript extends Feature<Exon> {
    public Transcript(IRegion region, String id) {
        super(region, id, null);
    }
}

@CompileStatic
class Exon extends Feature {
    int exonNumber
    Gene gene
    boolean coding
    public Exon(IRegion region, String id, int exonNumber) {
        super(region, id, null);
        this.exonNumber = exonNumber
    }
}
