package gngs.plot

import java.awt.geom.Rectangle2D
import java.awt.image.BufferedImage

import javax.imageio.ImageIO

import de.erichseifert.gral.plots.XYPlot

import org.junit.Test

/**
 * Covers sizing the image so that a legend placed outside the plot area is not
 * cut off, which is what the {@code marginRight} option and the legend width
 * measurement exist for.
 */
class LegendWidthTest {

    private static final String LONG_NAME = '20250207_NA12878_promethion_run3_barcode'

    private Plot plotWithLegend(String location, int series = 3) {
        Plot p = new Plot(title: 'Read length distribution', xBound: [0, 100], yBound: [0, 100])
        p.legendLocation = location
        series.times { int i ->
            p << new Lines(x: [0, 100], y: [0, 100], displayName: LONG_NAME + i)
        }
        return p
    }

    /**
     * The right hand edge the legend reaches, once the plot is laid out
     */
    private double legendMaxX(Plot p, int width, int height) {
        XYPlot xyPlot = p.toXYPlot(width, height)
        xyPlot.setBounds(0, 0, width, height)
        return xyPlot.getLegend().getBounds().getMaxX()
    }

    // ------------------------------------------------ case insensitivity --

    @Test
    void 'the legend width estimate does not care about the case of the location'() {
        int lower = plotWithLegend('east').estimateLegendWidth(1024)
        int upper = plotWithLegend('EAST').estimateLegendWidth(1024)
        int mixed = plotWithLegend('East').estimateLegendWidth(1024)

        assert lower > 0
        assert upper == lower
        assert mixed == lower

        assert plotWithLegend('NORTH_EAST').estimateLegendWidth(1024) > 0
        assert plotWithLegend('SOUTH_EAST').estimateLegendWidth(1024) > 0
    }

    @Test
    void 'a location with nothing to the east needs no estimate'() {
        assert plotWithLegend('WEST').estimateLegendWidth(1024) == 0
        assert plotWithLegend('NORTH').estimateLegendWidth(1024) == 0
    }

    @Test
    void 'no legend location at all is handled'() {
        Plot p = new Plot()
        p << new Lines(x: [0, 1], y: [0, 1])

        assert p.legendLocation == null
        assert p.estimateLegendWidth(1024) == 0
    }

    // ------------------------------------------------- measured overflow --

    @Test
    void 'the measured overflow is exactly what the legend needs'() {

        Plot p = plotWithLegend('EAST')

        XYPlot xyPlot = p.toXYPlot(1024, 800)
        xyPlot.setBounds(0, 0, 1024, 800)

        int overflow = p.legendOverflowWidth(xyPlot, 1024)
        double needed = xyPlot.getLegend().getBounds().getMaxX() - 1024

        assert overflow >= needed
        assert overflow - needed < 1.0d
    }

    @Test
    void 'a legend that does not overhang needs no margin'() {

        Plot p = plotWithLegend('WEST')

        XYPlot xyPlot = p.toXYPlot(1024, 800)
        xyPlot.setBounds(0, 0, 1024, 800)

        assert p.legendOverflowWidth(xyPlot, 1024) == 0
    }

    @Test
    void 'a plot with no legend needs no margin'() {

        Plot p = new Plot(xBound: [0, 10], yBound: [0, 10])
        p << new Lines(x: [0, 1], y: [0, 1])

        XYPlot xyPlot = p.toXYPlot(1024, 800)
        xyPlot.setBounds(0, 0, 1024, 800)

        assert !xyPlot.isLegendVisible()
        assert p.legendOverflowWidth(xyPlot, 1024) == 0
    }

    // ------------------------------------------------------ getImage options --

    @Test
    void 'getImage takes the same options as save'() {
        Plot p = plotWithLegend('EAST')

        assert p.getImage(width: 600, height: 400).height == 400
        assert p.getImage(width: 600, height: 400, marginRight: 300).width == 900
    }

    @Test
    void 'getImage and save default their size the same way'() {
        Plot p = plotWithLegend('EAST')

        File out = new File('test.legend.default.png')
        if(out.exists())
            out.delete()
        p.save(out.path)

        BufferedImage saved = ImageIO.read(out)
        BufferedImage shown = p.getImage()

        assert shown.width == saved.width
        assert shown.height == saved.height
        assert shown.height == 800
    }

    @Test
    void 'getImage honours initWidth and initHeight'() {
        Plot p = plotWithLegend('EAST')
        p.initWidth = 700
        p.initHeight = 500

        BufferedImage image = p.getImage()
        assert image.height == 500
        assert image.width >= 700
    }

    @Test
    void 'the width and height overload still works'() {
        Plot p = plotWithLegend('EAST')
        assert p.getImage(600, 400).height == 400
    }

    // --------------------------------------------- the legend fits by default --

    @Test
    void 'an east legend fits without asking for a margin'() {

        for(String location in ['EAST', 'east', 'NORTH_EAST', 'SOUTH_EAST']) {

            Plot p = plotWithLegend(location)

            BufferedImage image = p.getImage()
            double needed = legendMaxX(p, 1024, 800)

            assert image.width >= needed : "legend cut off in getImage with location $location"

            File out = new File('test.legend.fit.png')
            if(out.exists())
                out.delete()
            p.save(out.path)

            assert ImageIO.read(out).width >= needed : "legend cut off in save with location $location"
        }
    }

    @Test
    void 'an explicit marginRight is still honoured'() {
        Plot p = plotWithLegend('EAST')

        File out = new File('test.legend.margin.png')
        if(out.exists())
            out.delete()
        p.save(out.path, marginRight: 500)

        assert ImageIO.read(out).width == 1524
        assert p.getImage(marginRight: 500).width == 1524
    }

    @Test
    void 'a wide legend is accommodated however long the names are'() {

        Plot p = new Plot(title: 't', xBound: [0, 100], yBound: [0, 100])
        p.legendLocation = 'EAST'
        6.times { int i ->
            p << new Lines(x: [0, 100], y: [0, 100],
                displayName: "${LONG_NAME}_and_then_some_more_text_to_be_sure_$i")
        }

        BufferedImage image = p.getImage()
        assert image.width >= legendMaxX(p, 1024, 800)
    }
}
