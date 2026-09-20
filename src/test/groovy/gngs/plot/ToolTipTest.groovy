package gngs.plot

import java.awt.Color
import java.awt.geom.Point2D
import java.awt.geom.Rectangle2D
import java.awt.image.BufferedImage

import de.erichseifert.gral.graphics.DrawingContext
import de.erichseifert.gral.plots.XYPlot
import de.erichseifert.gral.plots.axes.Axis
import de.erichseifert.gral.plots.axes.AxisRenderer
import de.erichseifert.gral.util.PointND

import org.junit.Test

class ToolTipTest {

    /**
     * Build the plot, lay it out and draw it, returning the GRAL plot, the
     * rendered image and the tooltip layer.
     */
    private Map draw(Plot p, int width = 800, int height = 600) {

        XYPlot xyPlot = p.toXYPlot(width, height)

        BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB)
        DrawingContext context = PlotUtils.createDrawingContext(image)

        xyPlot.setBounds(0, 0, width, height)
        xyPlot.draw(context)

        return [
            xyPlot: xyPlot,
            image: image,
            layer: xyPlot.find { it instanceof ToolTipLayer }
        ]
    }

    /**
     * Independently compute where GRAL puts a data value, so that the layer's
     * own resolution can be checked against it
     */
    private Point2D gralPosition(XYPlot xyPlot, Number x, Number y) {

        Rectangle2D plotBounds = xyPlot.getPlotArea().getBounds()
        Axis axisX = xyPlot.getAxis(XYPlot.AXIS_X)
        Axis axisY = xyPlot.getAxis(XYPlot.AXIS_Y)
        AxisRenderer rendererX = xyPlot.getAxisRenderer(XYPlot.AXIS_X)
        AxisRenderer rendererY = xyPlot.getAxisRenderer(XYPlot.AXIS_Y)

        return new Point2D.Double(
            plotBounds.getMinX() + rendererX.getPosition(axisX, x, true, false).get(PointND.X),
            plotBounds.getMinY() + rendererY.getPosition(axisY, y, true, false).get(PointND.Y)
        )
    }

    private Plot examplePlot() {
        return new Plot(title: 'An example of showing a tooltip', xBound: [0, 5], yBound: [0, 10]) <<
            new Points(x: [1, 2, 3, 4], y: [5, 6, 7, 8], toolTip: [1, 2, 3, 4].collect { "X value is: $it" })
    }

    // ------------------------------------------------------------- placement --

    @Test
    void 'compass anchors map to the expected angles'() {
        assert ToolTipPlacement.from('E').angle == 0.0d
        assert ToolTipPlacement.from('NE').angle == 45.0d
        assert ToolTipPlacement.from('N').angle == 90.0d
        assert ToolTipPlacement.from('NW').angle == 135.0d
        assert ToolTipPlacement.from('W').angle == 180.0d
        assert ToolTipPlacement.from('SW').angle == 225.0d
        assert ToolTipPlacement.from('S').angle == 270.0d
        assert ToolTipPlacement.from('SE').angle == 315.0d
    }

    @Test
    void 'anchor names are case and separator insensitive'() {
        assert ToolTipPlacement.from('nw').angle == 135.0d
        assert ToolTipPlacement.from('NORTH_WEST').angle == 135.0d
        assert ToolTipPlacement.from('north-west').angle == 135.0d
        assert ToolTipPlacement.from('north west').angle == 135.0d
    }

    @Test
    void 'readable direction names are accepted'() {
        assert ToolTipPlacement.from('above').angle == 90.0d
        assert ToolTipPlacement.from('below').angle == 270.0d
        assert ToolTipPlacement.from('left').angle == 180.0d
        assert ToolTipPlacement.from('right').angle == 0.0d
    }

    @Test
    void 'null and true mean automatic placement'() {
        assert ToolTipPlacement.from(null).angle == null
        assert ToolTipPlacement.from(true).angle == null
    }

    @Test
    void 'a number is taken as an angle'() {
        assert ToolTipPlacement.from(30).angle == 30.0d
        assert ToolTipPlacement.from(30.5d).angle == 30.5d
    }

    @Test
    void 'a map may specify angle and distance'() {
        ToolTipPlacement p = ToolTipPlacement.from([angle: 20, distance: 60])
        assert p.angle == 20.0d
        assert p.distance == 60.0d

        ToolTipPlacement byAnchor = ToolTipPlacement.from([anchor: 'NW', distance: 40])
        assert byAnchor.angle == 135.0d
        assert byAnchor.distance == 40.0d

        ToolTipPlacement distanceOnly = ToolTipPlacement.from([distance: 40])
        assert distanceOnly.angle == null
        assert distanceOnly.distance == 40.0d
    }

    @Test
    void 'an unknown anchor is rejected with a helpful message'() {
        try {
            ToolTipPlacement.from('upwards')
            assert false : 'expected an exception'
        }
        catch(IllegalArgumentException e) {
            assert e.message.contains('upwards')
            assert e.message.contains('north')
        }
    }

    // ------------------------------------------------------------- selection --

    @Test
    void 'tooltips can be selected by a list of indices'() {
        Points points = new Points(x: [1, 2, 3], y: [1, 2, 3], toolTip: ['a', 'b', 'c'])
        points.showTooltips([0, 2])

        assert points.shownToolTips.keySet() as List == [0, 2]
        assert points.shownToolTips.values().every { it.angle == null }
    }

    @Test
    void 'tooltips can be selected by a single index'() {
        Points points = new Points(x: [1, 2, 3], y: [1, 2, 3], toolTip: ['a', 'b', 'c'])
        points.showTooltips(1)

        assert points.shownToolTips.keySet() as List == [1]
    }

    @Test
    void 'a map selects tooltips and places each one'() {
        Points points = new Points(x: [1, 2, 3], y: [1, 2, 3], toolTip: ['a', 'b', 'c'])
        points.showTooltips([0: 'NW', 2: [angle: 10, distance: 50]])

        assert points.shownToolTips[0].angle == 135.0d
        assert points.shownToolTips[2].angle == 10.0d
        assert points.shownToolTips[2].distance == 50.0d
    }

    @Test
    void 'a closure selects tooltips by data value'() {
        Points points = new Points(x: [1, 2, 3, 4], y: [10, 20, 30, 40], toolTip: ['a', 'b', 'c', 'd'])
        points.showTooltips { x, y, tip, i -> y > 25 }

        assert points.shownToolTips.keySet() as List == [2, 3]
    }

    @Test
    void 'a closure may return a placement rather than just true'() {
        Points points = new Points(x: [1, 2, 3], y: [1, 2, 3], toolTip: ['a', 'b', 'c'])
        points.showTooltips { x, y, tip, i -> i == 1 ? 'NW' : null }

        assert points.shownToolTips.keySet() as List == [1]
        assert points.shownToolTips[1].angle == 135.0d
    }

    @Test
    void 'a selection closure may declare fewer parameters'() {
        Points points = new Points(x: [1, 2, 3], y: [5, 6, 7], toolTip: ['a', 'b', 'c'])
        points.showTooltips { x -> x > 2 }

        assert points.shownToolTips.keySet() as List == [2]
    }

    @Test
    void 'showAllTooltips selects every point that has text'() {
        Points points = new Points(x: [1, 2, 3], y: [1, 2, 3], toolTip: ['a', null, 'c'])
        points.showAllTooltips()

        assert points.shownToolTips.keySet() as List == [0, 2]
    }

    // ---------------------------------------------------------- tooltip text --

    @Test
    void 'tooltips given as a list are resolved in order'() {
        Points points = new Points(x: [1, 2], y: [1, 2], toolTip: ['a', 'b'])

        assert points.resolveToolTips() == ['a', 'b']
    }

    @Test
    void 'GStrings in a tooltip list are resolved to strings'() {
        Points points = new Points(x: [1, 2], y: [1, 2], toolTip: [1, 2].collect { "X value is: $it" })

        assert points.resolveToolTips() == ['X value is: 1', 'X value is: 2']
    }

    @Test
    void 'a tooltip closure is called per point with as many arguments as it declares'() {
        Points points = new Points(x: [1, 2, 3], y: [10, 20, 30], displayName: 'Series')

        points.toolTip = { x -> "x=$x" }
        assert points.resolveToolTips() == ['x=1', 'x=2', 'x=3']

        points.toolTip = { x, y -> "$x/$y" }
        assert points.resolveToolTips() == ['1/10', '2/20', '3/30']

        points.toolTip = { x, y, i -> "$i:$x" }
        assert points.resolveToolTips() == ['0:1', '1:2', '2:3']
    }

    @Test
    void 'no tooltips resolves to an empty list'() {
        assert new Points(x: [1], y: [1]).resolveToolTips() == []
    }

    @Test
    void 'an unusable tooltip value is rejected'() {
        try {
            new Points(x: [1], y: [1], toolTip: 42).resolveToolTips()
            assert false : 'expected an exception'
        }
        catch(IllegalArgumentException e) {
            assert e.message.contains('Integer')
        }
    }

    // ------------------------------------------------------------- plot level --

    @Test
    void 'a series can be selected by display name'() {
        Plot p = new Plot() <<
            new Points(x: [1, 2], y: [1, 2], toolTip: ['a', 'b'], displayName: 'Bananas') <<
            new Points(x: [1, 2], y: [3, 4], toolTip: ['c', 'd'], displayName: 'Oranges')

        p.showTooltips('Oranges', [1])

        List<XYItem> xys = p.items.grep { it instanceof XYItem }
        assert xys[0].shownToolTips == null
        assert xys[1].shownToolTips.keySet() as List == [1]
    }

    @Test
    void 'a series can be selected by index'() {
        Plot p = new Plot() <<
            new Points(x: [1, 2], y: [1, 2], toolTip: ['a', 'b']) <<
            new Points(x: [1, 2], y: [3, 4], toolTip: ['c', 'd'])

        p.showTooltips(0, [1])

        List<XYItem> xys = p.items.grep { it instanceof XYItem }
        assert xys[0].shownToolTips.keySet() as List == [1]
        assert xys[1].shownToolTips == null
    }

    @Test
    void 'with no selector every series that has tooltips is annotated'() {
        Plot p = new Plot() <<
            new Points(x: [1, 2], y: [1, 2], toolTip: ['a', 'b']) <<
            new Points(x: [1, 2], y: [3, 4], toolTip: ['c', 'd'])

        p.showTooltips([0])

        assert p.items.grep { it instanceof XYItem }.every { it.shownToolTips.keySet() as List == [0] }
    }

    @Test
    void 'an unknown display name is rejected'() {
        Plot p = new Plot() << new Points(x: [1], y: [1], toolTip: ['a'], displayName: 'Bananas')

        try {
            p.showTooltips('Oranges', [0])
            assert false : 'expected an exception'
        }
        catch(IllegalArgumentException e) {
            assert e.message.contains('Oranges')
            assert e.message.contains('Bananas')
        }
    }

    @Test
    void 'showing tooltips when none are set is rejected'() {
        Plot p = new Plot() << new Points(x: [1], y: [1])

        try {
            p.showTooltips([0])
            assert false : 'expected an exception'
        }
        catch(IllegalStateException e) {
            assert e.message.contains('toolTip')
        }
    }

    @Test
    void 'an out of range tooltip index is reported when the plot is built'() {
        Plot p = new Plot() << new Points(x: [1, 2], y: [1, 2], toolTip: ['a', 'b'])
        p.showTooltips([5])

        try {
            p.toXYPlot(800, 600)
            assert false : 'expected an exception'
        }
        catch(IllegalArgumentException e) {
            assert e.message.contains('5')
        }
    }

    // --------------------------------------------------------------- geometry --

    @Test
    void 'no tooltip layer is added when nothing is selected'() {
        Map rendered = draw(examplePlot())

        assert rendered.layer == null
    }

    @Test
    void 'a tooltip anchors exactly where GRAL renders its data point'() {

        Plot p = examplePlot()
        p.showTooltips([1])

        Map rendered = draw(p)
        ToolTipLayer layer = rendered.layer

        assert layer != null
        assert layer.annotations.size() == 1

        // The annotated point is index 1, ie: x=2, y=6
        Point2D expected = gralPosition(rendered.xyPlot, 2, 6)
        Point2D actual = layer.annotations[0].anchor

        assert Math.abs(actual.x - expected.x) < 0.0001
        assert Math.abs(actual.y - expected.y) < 0.0001
    }

    @Test
    void 'the marker is painted at the anchor'() {

        Plot p = examplePlot()
        p.showTooltips([1])

        Map rendered = draw(p)
        Point2D anchor = rendered.layer.annotations[0].anchor

        Color painted = new Color(rendered.image.getRGB((int)Math.round(anchor.x), (int)Math.round(anchor.y)))
        Color expected = new ToolTipStyle().leaderColor

        assert Math.abs(painted.red - expected.red) < 20
        assert Math.abs(painted.green - expected.green) < 20
        assert Math.abs(painted.blue - expected.blue) < 20
    }

    @Test
    void 'an angle of 90 degrees puts the tooltip above the point and 270 below'() {

        Plot above = examplePlot()
        above.showTooltips([1: 'N'])
        ToolTipAnnotation annotationAbove = draw(above).layer.annotations[0]

        assert annotationAbove.bounds.centerY < annotationAbove.anchor.y
        assert Math.abs(annotationAbove.bounds.centerX - annotationAbove.anchor.x) < 1.0d

        Plot below = examplePlot()
        below.showTooltips([1: 'S'])
        ToolTipAnnotation annotationBelow = draw(below).layer.annotations[0]

        assert annotationBelow.bounds.centerY > annotationBelow.anchor.y
    }

    @Test
    void 'east and west place the tooltip to the right and left'() {

        Plot east = examplePlot()
        east.showTooltips([1: 'E'])
        ToolTipAnnotation annotationEast = draw(east).layer.annotations[0]

        assert annotationEast.bounds.minX > annotationEast.anchor.x

        Plot west = examplePlot()
        west.showTooltips([1: 'W'])
        ToolTipAnnotation annotationWest = draw(west).layer.annotations[0]

        assert annotationWest.bounds.maxX < annotationWest.anchor.x
    }

    @Test
    void 'an explicit distance is honoured'() {

        Plot near = examplePlot()
        near.showTooltips([1: [anchor: 'E', distance: 20]])
        ToolTipAnnotation annotationNear = draw(near).layer.annotations[0]

        Plot far = examplePlot()
        far.showTooltips([1: [anchor: 'E', distance: 80]])
        ToolTipAnnotation annotationFar = draw(far).layer.annotations[0]

        double gapNear = annotationNear.bounds.minX - annotationNear.anchor.x
        double gapFar = annotationFar.bounds.minX - annotationFar.anchor.x

        assert Math.abs(gapNear - 20.0d) < 1.0d
        assert Math.abs(gapFar - 80.0d) < 1.0d
    }

    @Test
    void 'an explicit anchor keeps its direction but steps out to avoid overlapping'() {

        // Several series annotated in the same direction from points at a similar
        // height: without stepping out, every box lands on top of the last one
        Plot p = new Plot(xBound: [0, 10], yBound: [0, 10])
        (0..3).each { int s ->
            p << new Points(
                x: [2 + s], y: [8],
                displayName: "Series $s",
                toolTip: ["a reasonably long tooltip for series $s"])
        }

        p.showTooltips([0: [anchor: 'south_east']])

        List<ToolTipAnnotation> annotations = draw(p).layer.annotations
        assert annotations.size() == 4

        for(int i = 0; i < annotations.size(); ++i) {
            // the requested direction is honoured for every one of them
            assert annotations[i].bounds.centerX > annotations[i].anchor.x
            assert annotations[i].bounds.centerY > annotations[i].anchor.y

            for(int j = i + 1; j < annotations.size(); ++j) {
                Rectangle2D overlap = annotations[i].bounds.createIntersection(annotations[j].bounds)
                assert overlap.width <= 0 || overlap.height <= 0 :
                    "tooltip $i overlaps tooltip $j"
            }
        }
    }

    @Test
    void 'an explicit anchor is left exactly where asked when nothing is in the way'() {

        Plot p = examplePlot()
        p.showTooltips([1: [anchor: 'E', distance: 30]])

        ToolTipAnnotation annotation = draw(p).layer.annotations[0]

        assert Math.abs((annotation.bounds.minX - annotation.anchor.x) - 30.0d) < 1.0d
    }

    @Test
    void 'automatic placement keeps tooltips on adjacent points apart'() {

        Plot p = new Plot(xBound: [0, 10], yBound: [0, 10]) <<
            new Points(x: [4, 5], y: [5, 5], toolTip: ['first point here', 'second point here'])

        p.showAllTooltips()

        List<ToolTipAnnotation> annotations = draw(p).layer.annotations
        assert annotations.size() == 2

        Rectangle2D overlap = annotations[0].bounds.createIntersection(annotations[1].bounds)
        assert overlap.width <= 0 || overlap.height <= 0
    }

    @Test
    void 'automatic placement stays within the plot area'() {

        Plot p = new Plot(xBound: [0, 10], yBound: [0, 10]) <<
            new Points(x: [5], y: [5], toolTip: ['a tooltip in the middle'])

        p.showAllTooltips()

        Map rendered = draw(p)
        Rectangle2D plotBounds = rendered.xyPlot.getPlotArea().getBounds()

        assert plotBounds.contains(rendered.layer.annotations[0].bounds)
    }

    // ---------------------------------------------------------------- styling --

    @Test
    void 'style options passed to showTooltips apply to that series only'() {

        Plot p = new Plot() <<
            new Points(x: [1, 2], y: [1, 2], toolTip: ['a', 'b'], displayName: 'Bananas') <<
            new Points(x: [1, 2], y: [3, 4], toolTip: ['c', 'd'], displayName: 'Oranges')

        p.showTooltips('Oranges', [0], distance: 55, leader: false)

        List<XYItem> xys = p.items.grep { it instanceof XYItem }
        assert xys[0].toolTipStyle == null
        assert xys[1].toolTipStyle.distance == 55.0d
        assert !xys[1].toolTipStyle.leader
    }

    @Test
    void 'an unknown style option is rejected'() {
        Points points = new Points(x: [1], y: [1], toolTip: ['a'])

        try {
            points.showTooltips([0], colour: Color.red)
            assert false : 'expected an exception'
        }
        catch(MissingPropertyException e) {
            assert e.message.contains('colour')
        }
    }

    @Test
    void 'maxWidth wraps text into a narrower and taller tooltip'() {

        String text = 'a fairly long tooltip that will certainly need to be wrapped somewhere'

        Plot wide = new Plot(xBound: [0, 10], yBound: [0, 10]) << new Points(x: [5], y: [5], toolTip: [text])
        wide.showAllTooltips()
        Rectangle2D wideBounds = draw(wide).layer.annotations[0].bounds

        Plot narrow = new Plot(xBound: [0, 10], yBound: [0, 10]) << new Points(x: [5], y: [5], toolTip: [text])
        narrow.showAllTooltips(maxWidth: 120.0d)
        Rectangle2D narrowBounds = draw(narrow).layer.annotations[0].bounds

        assert narrowBounds.width < wideBounds.width
        assert narrowBounds.width <= 120.0d
        assert narrowBounds.height > wideBounds.height
    }

    @Test
    void 'markup is rendered rather than shown literally'() {

        Plot p = new Plot(xBound: [0, 10], yBound: [0, 10]) <<
            new Points(x: [5], y: [5], toolTip: ['<b>bold</b><br>second line'])

        p.showAllTooltips()
        Rectangle2D bounds = draw(p).layer.annotations[0].bounds

        Plot plain = new Plot(xBound: [0, 10], yBound: [0, 10]) <<
            new Points(x: [5], y: [5], toolTip: ['bold'])
        plain.showAllTooltips()
        Rectangle2D plainBounds = draw(plain).layer.annotations[0].bounds

        // Two lines of text, so taller than the equivalent single line
        assert bounds.height > plainBounds.height
    }

    // ------------------------------------------------------------ integration --

    @Test
    void 'the documented example saves a plot with tooltips'() {

        Plot p = new Plot(title: 'An example of showing a tooltip') <<
            new Points(x: [1, 2, 3, 4], y: [5, 6, 7, 8], toolTip: [1, 2, 3, 4].collect { "X value is: $it" })

        p.showTooltips([0, 2])

        File out = new File('test.tooltips.png')
        if(out.exists())
            out.delete()

        p.save(out.path)

        assert out.exists()
        assert out.length() > 0

        BufferedImage image = javax.imageio.ImageIO.read(out)
        assert image.width == 1024
        assert image.height == 800
    }

    @Test
    void 'tooltips render on a multi series plot with lines and markup'() {

        Plot p = new Plot(
            title: 'Tooltips with markup and automatic placement',
            xLabel: 'Position',
            yLabel: 'Depth',
            xBound: [0, 11],
            yBound: [0, 11]
        ) <<
            new Lines(
                x: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                y: [5, 6.2, 7, 6.4, 8.1, 8.6, 7.2, 9, 9.4, 8.2],
                displayName: 'Bananas',
                toolTip: (1..10).collect { "<b>point $it</b><br>depth: <i>${it * 2}x</i>" }
            ) <<
            new Lines(
                x: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                y: [2, 3.1, 2.6, 4, 3.4, 5.2, 4.6, 5.9, 5.1, 6.4],
                displayName: 'Oranges',
                toolTip: (1..10).collect { "sample &amp; $it" }
            )

        p.showTooltips('Bananas', [1: true, 5: true, 8: 'SW'])
        p.showTooltips('Oranges', [3, 4, 9])

        File out = new File('test.tooltips.multi.png')
        if(out.exists())
            out.delete()

        p.save(out.path)
        assert out.exists()
        assert out.length() > 0
    }

    @Test
    void 'tooltips transfer from a BeakerX points item'() {

        def bxPoints = new com.twosigma.beakerx.chart.xychart.plotitem.Points(
            x: [1, 2, 3], y: [4, 5, 6], displayName: 'From BeakerX', toolTip: ['a', 'b', 'c'])

        // Plot.from walks bxPlot.graphics, which is what this stands in for
        Points item = new Points()
        item.toolTip = Plot.beakerXToolTips(bxPoints)

        assert item.resolveToolTips() == ['a', 'b', 'c']
    }
}
