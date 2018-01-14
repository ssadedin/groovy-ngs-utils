package gngs;

import java.util.Iterator;
import java.util.List;
import java.util.Map;

import groovy.lang.IntRange;

public class RegionIterator implements Iterator<Region> {
        
    Iterator<Map.Entry<String,List<IntRange>>> allIterator = null;

    // Iterator chrIterator = allIterator.hasNext() ? allIterator.next().value.iterator() : null
    Iterator<IntRange> chrIterator =  null;

    Regions regions;

    RegionIterator(Regions regions) {
        this.regions = regions;
        Map<String, List<IntRange>> allRanges = this.regions.getAllRanges();
        allIterator = allRanges.entrySet().iterator();
    }

    Map.Entry<String, List<IntRange>> currentChr;

    String chr;

    public boolean hasNext() {
        if(chr == null) {
            for(List<IntRange> r : this.regions.getAllRanges().values()) {
                if(!r.isEmpty())
                    return true;
            }
            return false;
        }
        else
            return chrIterator.hasNext();
    }

    public Region next() {
        if(chr == null)
            nextChr();

        IntRange nextRange = chrIterator.next();
        String currentChr = chr;
        if(!chrIterator.hasNext())
            nextChr();

        if(nextRange instanceof GRange) {
            GRange grange = (GRange) nextRange;
            if(grange.getExtra() instanceof Region) {
                return (Region)grange.getExtra();
            }
        }
        return new Region(currentChr,nextRange);
    }

    void nextChr() {
        while(allIterator.hasNext()) {
            currentChr = allIterator.next();
            if(currentChr.getValue() != null) {
                chrIterator = currentChr.getValue().iterator();
                chr = currentChr.getKey();
                break;
            }
        }
    }

    public void remove() {
        throw new UnsupportedOperationException();
    }
}