package gngs.gencode

import com.google.common.collect.ImmutableList
import gngs.GRange
import gngs.IRegion
import groovy.transform.CompileStatic

@CompileStatic
class Feature implements IRegion {
    String chr
    IntRange range
    String id
    Feature parent
    List<Feature> children
    
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
        if(this.children.isEmpty())
            return Collections.emptyList()
        return this.children
    }
    
    String toString() {
        "${this.class.name.replace('gngs.gencode.','')} at $chr:$range.from-$range.to"
    }
}

@CompileStatic
class Gene extends Feature {
    String symbol
    Transcript transcript
    public Gene(IRegion region, String id, String symbol) {
        super(region, id, null);
        this.symbol = symbol
    }
}

@CompileStatic
class Transcript extends Feature {
    public Transcript(IRegion region, String id) {
        super(region, id, null);
    }
}

@CompileStatic
class Exon extends Feature {
    int exonNumber
    Gene gene
    public Exon(IRegion region, String id) {
        super(region, id, null);
    }
}
