import static org.junit.Assert.*;

import org.junit.Test;

class TargetedCNVAnnotatorTest {

    @Test
    public void testMaxRange() {
        
        RangedData dgv = new RangedData()
        
        Regions target = createTarget([
         100..120,
         140..160,
         200..220
        ])
        
        TargetedCNVAnnotator a = new TargetedCNVAnnotator(target, dgv)
        Variant v = del(140,20)
        def max = a.computeMaxRange(v)
        assert max.from == 120
        assert max.to == 200
    }
    
    @Test
    public void testCompatibleRanges() {
        
        RangedData d = dgv([125..165])
        
        Regions target = createTarget([
         100..120,
         140..160,
         200..220
        ])
        
        TargetedCNVAnnotator a = new TargetedCNVAnnotator(target, d)
        Variant v = del(140,20)
        
        def ranges = a.findCompatibleRanges(v)
        assert ranges.size() == 1
        
        a = new TargetedCNVAnnotator(target, this.dgv([115..165]))
        ranges = a.findCompatibleRanges(v)
        assert ranges.size() == 0
        
        a = new TargetedCNVAnnotator(target, this.dgv([125..205]))
        ranges = a.findCompatibleRanges(v)
        assert ranges.size() == 0
        
        a = new TargetedCNVAnnotator(target, this.dgv([125..155]))
        ranges = a.findCompatibleRanges(v)
        assert ranges.size() == 0
    }
    
    RangedData dgv(List<IntRange> ranges) {
        RangedData dgv = new RangedData()
        ranges.each { dgv.addRegion("chr1", it.from, it.to) }
        return dgv
    }
    
    Regions createTarget(List<IntRange> target) {
        Regions result = new Regions()
        target.each { result.addRegion("chr1", it.from, it.to) }
        return result
    }
    
    Variant del(int pos, int size) {
        Variant v = new Variant()
        v.chr = "chr1"
        v.pos = pos
        v.alt = "<DEL>"
        v.update { v.info.SVLEN=-size; v.info.END=pos+size }
        return v
    }
}