package gngs.plot

import java.awt.Color

import com.twosigma.beakerx.chart.xychart.plotitem.StrokeType

import gngs.plot.bx.Density

import org.junit.Test

/**
 * Covers adding BeakerX chart items to a gngs {@link Plot} with the left shift
 * operator, which is the path taken by notebook code that is moved into a script.
 */
class BeakerXItemsTest {

    private static com.twosigma.beakerx.chart.xychart.plotitem.Line bxLine(Map attributes = [:]) {
        return new com.twosigma.beakerx.chart.xychart.plotitem.Line([x: [1, 2, 3], y: [4, 5, 6]] + attributes)
    }

    // --------------------------------------------------------- new overloads --

    @Test
    void 'a BeakerX Points can be added to a plot'() {
        Plot p = new Plot()
        p << new com.twosigma.beakerx.chart.xychart.plotitem.Points(
            x: [1, 2, 3], y: [4, 5, 6], displayName: 'Bananas')

        assert p.items.size() == 1
        assert p.items[0] instanceof Points
        assert p.items[0].x == [1, 2, 3]
        assert p.items[0].y == [4, 5, 6]
        assert p.items[0].displayName == 'Bananas'
    }

    @Test
    void 'a BeakerX Bars can be added to a plot'() {
        Plot p = new Plot()
        p << new com.twosigma.beakerx.chart.xychart.plotitem.Bars(
            x: [1, 2, 3], y: [4, 5, 6], displayName: 'Oranges', width: 3)

        assert p.items.size() == 1
        assert p.items[0] instanceof Bars
        assert p.items[0].displayName == 'Oranges'
        assert p.items[0].width == 3.0d
    }

    @Test
    void 'a plot of BeakerX Points and Bars renders'() {
        Plot p = new Plot(title: 'from beakerx')
        p << new com.twosigma.beakerx.chart.xychart.plotitem.Bars(x: [1, 2, 3], y: [4, 5, 6], displayName: 'b')

        File out = new File('test.bx.bars.png')
        if(out.exists())
            out.delete()

        p.save(out.path)
        assert out.exists() && out.length() > 0
    }

    // ------------------------------------------------------------ the fix --

    @Test
    void 'the style of a BeakerX Line is carried across'() {
        Plot p = new Plot()
        p << bxLine(style: StrokeType.DASH)

        // gngs holds the style as a string, which is what the renderer matches on
        assert p.items[0].style == 'DASH'
    }

    @Test
    void 'every stroke type is carried across by name'() {
        for(StrokeType type in StrokeType.values()) {
            Plot p = new Plot()
            p << bxLine(style: type)
            assert p.items[0].style == type.name()
        }
    }

    @Test
    void 'a dashed BeakerX line renders as a dashed line'() {
        Plot p = new Plot(xBound: [0, 4], yBound: [0, 8])
        p << bxLine(style: StrokeType.DASH, width: 2)

        assert p.items[0].style == 'DASH'

        File out = new File('test.bx.dashed.png')
        if(out.exists())
            out.delete()
        p.save(out.path)
        assert out.exists() && out.length() > 0
    }

    // ------------------------------------------------------ no regressions --

    @Test
    void 'the attributes previously copied by hand are still copied'() {
        Plot p = new Plot()
        p << bxLine(displayName: 'Bananas', width: 4)

        assert p.items[0] instanceof Line
        assert p.items[0].x == [1, 2, 3]
        assert p.items[0].y == [4, 5, 6]
        assert p.items[0].displayName == 'Bananas'
        assert p.items[0].width == 4.0d
    }

    @Test
    void 'a BeakerX colour becomes an awt colour'() {
        Plot p = new Plot()
        p << bxLine(color: com.twosigma.beakerx.chart.Color.RED)

        assert p.items[0].color instanceof Color
        assert p.items[0].color.red == 255
        assert p.items[0].color.green == 0
        assert p.items[0].color.blue == 0
    }

    @Test
    void 'an unset colour is left for the palette to fill in'() {
        Plot p = new Plot()
        p << bxLine()

        assert p.items[0].color == null

        // which means the palette still assigns by series position
        assert p.convertColor(p.items[0].color, 1) == new DefaultPalette().colors[1]
    }

    @Test
    void 'a BeakerX Area still converts'() {
        Plot p = new Plot()
        p << new com.twosigma.beakerx.chart.xychart.plotitem.Area(
            x: [1, 2, 3], y: [4, 5, 6], displayName: 'Area')

        assert p.items[0] instanceof Area
        assert p.items[0].x == [1, 2, 3]
        assert p.items[0].displayName == 'Area'
    }

    @Test
    void 'a density area still converts and renders'() {
        // Density subclasses the BeakerX types and computes x and y from data,
        // so it goes through the same conversion and must not be disturbed by it
        Random r = new Random(11)
        List values = (1..500).collect { r.nextGaussian() }

        Plot p = new Plot()
        Density.Area area = new Density.Area(
            data: values, displayName: 'Density', color: com.twosigma.beakerx.chart.Color.blue)
        p << area

        assert p.items.size() == 1
        assert p.items[0] instanceof Area
        assert p.items[0].displayName == 'Density'
        assert p.items[0].x.size() > 0
        assert p.items[0].x.size() == p.items[0].y.size()
        assert p.items[0].color instanceof Color

        File out = new File('test.bx.density.png')
        if(out.exists())
            out.delete()
        p.save(out.path)
        assert out.exists() && out.length() > 0
    }

    @Test
    void 'properties the gngs item does not have are ignored'() {
        // BeakerX carries plenty that has no gngs equivalent, eg: lodFilter,
        // plotType, outlineColor. None of it should cause trouble.
        Plot p = new Plot()
        p << new com.twosigma.beakerx.chart.xychart.plotitem.Points(
            x: [1, 2], y: [3, 4], size: 12, shape: com.twosigma.beakerx.chart.xychart.plotitem.ShapeType.CIRCLE)

        assert p.items.size() == 1
        assert p.items[0].x == [1, 2]
    }

    // ----------------------------------------- left shift versus Plot.from --

    @Test
    void 'left shift and Plot_from share one conversion'() {

        // Plot.from cannot be exercised here because a BeakerX Plot needs a live
        // kernel to construct. Both paths now call the same copier though, so
        // driving it directly covers the conversion they share, including the
        // Lines class that Plot.from creates where left shift creates a Line.
        def source = bxLine(displayName: 'Bananas', width: 2, style: StrokeType.DOT,
                            color: com.twosigma.beakerx.chart.Color.GREEN)

        Plot shifted = new Plot()
        shifted << source
        XYItem viaShift = (XYItem)shifted.items[0]

        Lines viaFrom = new Lines()
        new Plot().copyBeakerXProperties(source, viaFrom, 0)

        assert viaShift.x == viaFrom.x
        assert viaShift.y == viaFrom.y
        assert viaShift.displayName == viaFrom.displayName
        assert viaShift.width == viaFrom.width
        assert viaShift.style == viaFrom.style
        assert viaShift.color == viaFrom.color

        assert viaFrom.style == 'DOT'
        assert viaFrom.color instanceof Color
    }

    @Test
    void 'tooltips are still not carried across by left shift'() {
        // Documents a known gap rather than endorsing it: BeakerX exposes
        // toolTips while gngs calls it toolTip, so the names do not match
        Plot p = new Plot()
        p << bxLine(toolTip: ['a', 'b', 'c'])

        assert p.items[0].toolTip == null
    }
}
