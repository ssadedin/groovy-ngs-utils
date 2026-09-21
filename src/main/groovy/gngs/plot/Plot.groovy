package gngs.plot

import java.awt.BasicStroke
import java.awt.Color
import java.awt.Graphics2D
import java.awt.Paint
import java.awt.RenderingHints
import java.awt.Stroke
import java.awt.geom.Rectangle2D
import java.awt.image.BufferedImage

import javax.imageio.ImageIO
import javax.imageio.ImageWriter
import javax.imageio.stream.ImageOutputStream

import com.twosigma.beakerx.chart.xychart.plotitem.XYGraphics

import de.erichseifert.gral.data.Column
import de.erichseifert.gral.data.DataSource
import de.erichseifert.gral.data.DataTable
import de.erichseifert.gral.data.EnumeratedData
import de.erichseifert.gral.data.statistics.Statistics
import de.erichseifert.gral.data.statistics.Histogram2D
import de.erichseifert.gral.graphics.Drawable
import de.erichseifert.gral.graphics.DrawingContext
import de.erichseifert.gral.graphics.DrawingContext.Quality
import de.erichseifert.gral.graphics.DrawingContext.Target
import de.erichseifert.gral.graphics.Insets2D
import de.erichseifert.gral.graphics.Label
import de.erichseifert.gral.graphics.Location
import de.erichseifert.gral.graphics.Orientation
import de.erichseifert.gral.io.plots.DrawableWriter
import de.erichseifert.gral.io.plots.DrawableWriterFactory
import de.erichseifert.gral.plots.BarPlot
import de.erichseifert.gral.plots.BarPlot.BarPlotLegend
import de.erichseifert.gral.plots.BarPlot.BarRenderer
import de.erichseifert.gral.plots.XYPlot
import de.erichseifert.gral.plots.XYPlot.XYLegend
import de.erichseifert.gral.plots.areas.AreaRenderer
import de.erichseifert.gral.plots.areas.DefaultAreaRenderer2D
import de.erichseifert.gral.plots.axes.Axis
import de.erichseifert.gral.plots.axes.AxisRenderer
import de.erichseifert.gral.plots.legends.Legend
import de.erichseifert.gral.plots.legends.SeriesLegend
import de.erichseifert.gral.plots.lines.DefaultLineRenderer2D
import de.erichseifert.gral.plots.lines.DiscreteLineRenderer2D
import de.erichseifert.gral.plots.lines.LineRenderer
import de.erichseifert.gral.plots.lines.SmoothLineRenderer2D
import de.erichseifert.gral.plots.points.DefaultPointRenderer2D
import de.erichseifert.gral.plots.points.PointRenderer
import de.erichseifert.gral.util.GraphicsUtils
import graxxia.Stats
import groovy.transform.CompileStatic


class Palette {
    Color [] colors 
    
}

class DefaultPalette extends Palette {
    DefaultPalette() {
       colors = [new Color(39,119,180), new Color(255, 127,14), new Color(44,160,44), Color.red, new Color(148,103,189), new Color(100,0,100), Color.cyan, Color.pink, Color.yellow, Color.magenta ] 
    }
}

class PlotItem {
    String displayName = null
}

class ConstantLine {
    Double x
    Double y
    Double width
    Object color
    String style
    String displayName
    
    DataTable toTable(double minX, double maxX, double minY, double maxY) {
        Column xColumn
        Column yColumn

        if(this.x != null) {
            xColumn = new Column(Double, [this.x, this.x])
            yColumn = new Column(Double, [minY, maxY])
        }
        else
        if(this.y != null) {
            xColumn = new Column(Double, [minX, maxX])
            yColumn = new Column(Double, [this.y, this.y])
        }

        return new DataTable(xColumn, yColumn)
    }
}

class Text {
    String text
    double x
    double y
    Color color
}

class XYItem extends PlotItem {
    Iterable<Object> x
    Iterable<Number> y
    
    Object color

    /**
     * Tooltip text for each data point, in the same order as {@link #x} and {@link #y}.
     * <p>
     * Accepts either a list of values, or - as BeakerX does - a closure which is
     * called for each point to build its text. The closure may declare any of
     * {@code (x, y, index, base, displayName)}, and only as many arguments as it
     * declares are passed.
     * <p>
     * Assigning tooltips does not by itself display them: see
     * {@link #showTooltips(Map,Object)}.
     */
    Object toolTip

    /**
     * Which tooltips to render, and where to put each one, keyed by the index of
     * the data point within this series.
     * <p>
     * This is an extension beyond the BeakerX interface, which has no equivalent
     * because it renders tooltips interactively.
     */
    Map<Integer, ToolTipPlacement> shownToolTips = null

    /**
     * Style used for this series' tooltips. If null, the style of the enclosing
     * {@link Plot} is used.
     */
    ToolTipStyle toolTipStyle = null

    double maxX = Double.MIN_VALUE
    double maxY = Double.MIN_VALUE
    
    double minX = Double.MAX_VALUE
    double minY = Double.MAX_VALUE
    
     DataTable asType(Class clazz) {
        if(clazz == DataTable) {
            return toTable()
        }
    }
    
    DataTable toTable() {
        List xList = x as List
        List yList = y as List

        Column xColumn = new Column(Double, xList)
        Column yColumn = new Column(Double, yList)
        
        List<Column> cols = [xColumn, yColumn]
        
        minX = x.min()
        minY = y.min()
        maxX = x.max()
        maxY = y.max()

        DataTable dt =  createTable(cols)
        return dt ; 
    }
    
    DataTable createTable(List<Column> columns) {
        return new DataTable(*columns)
    }

    /**
     * Resolve {@link #toolTip} to one string per data point.
     *
     * @return list of tooltip text, which may be shorter than the data if no
     *         tooltips were set, or contain nulls for points with no tooltip
     */
    List<String> resolveToolTips() {

        if(toolTip == null)
            return []

        if(toolTip instanceof Closure) {

            Closure builder = (Closure)toolTip
            List xList = x as List
            List yList = y as List
            int argCount = Math.min(builder.maximumNumberOfParameters, 5)

            return (0..<xList.size()).collect { int i ->
                List args = [xList[i], yList[i], i, null, displayName]
                Object result = builder.call(*args[0..<argCount])
                return result?.toString()
            }
        }

        if(toolTip instanceof Iterable)
            return ((Iterable)toolTip).collect { it?.toString() }

        throw new IllegalArgumentException(
            'toolTip should be a list of values or a closure, but was: ' + toolTip.class.name)
    }

    /**
     * Display the tooltips at the given data points.
     * <p>
     * The points to annotate may be given as:
     * <ul>
     *   <li>a list of indices, eg: {@code [1,3]}, placed automatically</li>
     *   <li>a map of index to placement, eg: {@code [1: 'NW', 3: [angle:20, distance:60]]}</li>
     *   <li>a single index</li>
     *   <li>a closure, called with any of {@code (x, y, toolTip, index)}, returning
     *       a falsy value to skip the point, {@code true} to place it automatically,
     *       or a placement such as {@code 'NW'}</li>
     * </ul>
     * Placements are described in {@link ToolTipPlacement#from(Object)}.
     *
     * @param style optional named arguments overriding {@link ToolTipStyle} properties
     * @param which which tooltips to show
     */
    XYItem showTooltips(Map style = null, Object which) {

        if(style)
            this.toolTipStyle = (this.toolTipStyle ?: new ToolTipStyle()).copy(style)

        if(this.shownToolTips == null)
            this.shownToolTips = new LinkedHashMap<Integer, ToolTipPlacement>()

        this.shownToolTips.putAll(selectToolTips(which))

        return this
    }

    /**
     * Display every tooltip that has been set on this series
     *
     * @param style optional named arguments overriding {@link ToolTipStyle} properties
     */
    XYItem showAllTooltips(Map style = null) {
        List<String> tips = resolveToolTips()
        return showTooltips(style?:[:], (0..<tips.size()).grep { int i -> tips[i] } )
    }

    /**
     * Interpret a tooltip selection into indices and their placements
     */
    private Map<Integer, ToolTipPlacement> selectToolTips(Object which) {

        Map<Integer, ToolTipPlacement> result = new LinkedHashMap<Integer, ToolTipPlacement>()

        if(which instanceof Closure) {

            Closure filter = (Closure)which
            List xList = x as List
            List yList = y as List
            List<String> tips = resolveToolTips()
            int argCount = Math.min(filter.maximumNumberOfParameters, 4)

            for(int i = 0; i < xList.size(); ++i) {
                List args = [xList[i], yList[i], i < tips.size() ? tips[i] : null, i]
                Object outcome = filter.call(*args[0..<argCount])
                if(outcome)
                    result[i] = ToolTipPlacement.from(outcome instanceof Boolean ? null : outcome)
            }
        }
        else
        if(which instanceof Map) {
            ((Map)which).each { Object index, Object placement ->
                result[((Number)index).intValue()] = ToolTipPlacement.from(placement)
            }
        }
        else
        if(which instanceof Number) {
            result[((Number)which).intValue()] = new ToolTipPlacement()
        }
        else
        if(which instanceof Iterable) {
            ((Iterable)which).each { Object index ->
                result[((Number)index).intValue()] = new ToolTipPlacement()
            }
        }
        else {
            throw new IllegalArgumentException(
                'Tooltips to show should be given as a list of indices, a map of index to placement, ' +
                'or a closure, but was: ' + which?.class?.name)
        }

        return result
    }
}

class Lines extends XYItem {
    Double width
    String style
}

class Line extends XYItem {
    Double width
    String style
}

class Area extends XYItem {
}

class Points extends XYItem {
}



class Bars extends XYItem {
    Double width
    List labels
    
    DataTable createTable(List<Column> columns) {
        if(labels) {
            columns << new Column(String, labels*.toString())
        }
        if(width == null && x.size()>1) {
            width = (0.84) * Math.abs(x[1] - x[0])
        }
        return super.createTable(columns)
    }
}

class Histogram {
    
    String title
    
    /**
     * Accept either Iterable<Double> or Iterable<Iterable<Dobule>>
     */
    def data
    
    int binCount = 10
    
    Palette palette = new DefaultPalette()
    
    String xLabel
    
    String yLabel
    
    Double rangeMax = Double.POSITIVE_INFINITY

    Double rangeMin = Double.NEGATIVE_INFINITY
    
    Iterable<String> names
    
    BarPlot createPlot() {

        Iterable<Iterable<Double>> values
        if(data[0] instanceof Iterable) {
            values = data
        }
        else {
            values = [data]
        }
        
        List breaks = calculateBreaks(values)
        
        double width = breaks[1] - breaks[0]
        
        List nameList = names ? names as List : [null] * values.size()
         
        List<DataSource> histogram2ds = [values,nameList].transpose().collect { valueAndName ->
            createHistogramDataSource(valueAndName[1], width, valueAndName[0], breaks)
        }
        
        // Create new bar plot
        BarPlot plot = new BarPlot(*histogram2ds);
        plot.setBarWidth(width*0.84)
        
        plot.setBackground(Color.white)
        plot.getTitle().setText(title)        


        plot.getAxisRenderer(XYPlot.AXIS_X).with { 
            if(xLabel)
                label.text = xLabel
        }
        
        plot.getAxisRenderer(XYPlot.AXIS_Y).with { 
            if(yLabel)
                label.text = yLabel
            intersection = -Double.MAX_VALUE // left align the axis
        }

        if(rangeMin != Double.NEGATIVE_INFINITY || rangeMax != Double.POSITIVE_INFINITY) {
            plot.getAxis(XYPlot.AXIS_X).with {
                if(rangeMin != Double.NEGATIVE_INFINITY)
                    min = rangeMin - (breaks[1] - breaks[0])
                if(rangeMax != Double.POSITIVE_INFINITY)
                    max = rangeMax + (breaks[1] - breaks[0])
            }
        }

        plot.getAxis(XYPlot.AXIS_Y).with {
            double yMax = 0
            for(DataSource ds in histogram2ds) {
                for(int row = 0; row < ds.getRowCount(); ++row) {
                    double val = ((Number)ds.get(1, row)).doubleValue()
                    if(val > yMax)
                        yMax = val
                }
            }
            max = PlotUtils.roundUpToOOM(yMax)
        }
        
        histogram2ds.eachWithIndex { h2d, i ->
            PointRenderer barRenderer = plot.getPointRenderers(h2d).get(0);
            barRenderer.setColor(GraphicsUtils.deriveWithAlpha(palette.colors[i+1], 128));
        }
       
        Insets2D.Double insets = new Insets2D.Double(40.0, 80.0, 80.0, 80.0)
        plot.setInsets(insets);
        
        return plot
    }
    
    BufferedImage getImage() {
        getImage(800,600)
    }

    BufferedImage getImage(int width, int height) {
        BarPlot plot = createPlot()
        BufferedImage bImage = new BufferedImage(1024, 800, BufferedImage.TYPE_INT_ARGB);
        DrawingContext context = PlotUtils.createDrawingContext(bImage)
        plot.setBounds(0, 0, width, height);
        plot.draw(context)
        drawHistogramLegend(plot, context)
        return bImage
     }
    
    void save(final String fileName) {
        
        BarPlot plot = createPlot()

        // Unfortunately this results in a poor legend intended for something more
        // like a categoryh plot because it contains every value
        // plot.setLegendVisible(true)
        
        new File(fileName).withOutputStream { w ->
            DrawableWriter wr = DrawableWriterFactory.getInstance().get("image/png");
            PlotUtils.write(plot, w, 0,0, 1024, 800, 0) { DrawingContext ctx ->
                drawHistogramLegend(plot, ctx)
            }
        }        
    }
    
    void drawHistogramLegend(BarPlot plot, DrawingContext ctx) {
        Graphics2D g = ctx.getGraphics()
        if(names) {
            names.eachWithIndex { name, index ->
                g.setFont(plot.getFont().deriveFont(80))
                g.setColor(GraphicsUtils.deriveWithAlpha(palette.colors[index+1], 128))
                g.drawString(name, 100, 120 + index * 36)
            }
        }
    }
    
    DataSource createHistogramDataSource(String name, Double width, Iterable<Double> data, List breaks) {

        DataTable dt = new DataTable(1, Double)
        
        for(double d in data) {
            if(d >= rangeMin && d <= rangeMax)
                dt.add(d)
        }
        
        assert dt.columnCount == 1
        
        if(name != null)
            dt.setName(name)
       
        // Create histogram from data
        Histogram2D hist = new Histogram2D(dt, Orientation.VERTICAL, [breaks as Double[]] as Double[][]);
        
        // Create a second dimension (x axis) for plotting
        DataSource histogram2d = new EnumeratedData(hist, breaks.min()+width/2, width);    
        
        return histogram2d
    }

    static Histogram from(gngs.plot.Histogram other) {
        return other
    }

    /**
     * Convert a BeakerX Histogram to a gngs Histogram
     */
    static Histogram from(com.twosigma.beakerx.chart.histogram.Histogram bxHist) {
        Histogram h = new Histogram()
        h.title = bxHist.title ?: ''

        h.data = bxHist.getData()

        if(bxHist.getBinCount())
            h.binCount = bxHist.getBinCount()

        if(bxHist.getXLabel())
            h.xLabel = bxHist.getXLabel()

        if(bxHist.getYLabel())
            h.yLabel = bxHist.getYLabel()

        if(bxHist.getNames() && !bxHist.getNames().isEmpty())
            h.names = bxHist.getNames()

        if(bxHist.getRangeMin() != null)
            h.rangeMin = bxHist.getRangeMin().toDouble()

        if(bxHist.getRangeMax() != null)
            h.rangeMax = bxHist.getRangeMax().toDouble()

        if(bxHist.getColor() && !bxHist.getColor().isEmpty()) {
            List<Color> colors = bxHist.getColor().collect { c ->
                new Color(c.getRed(), c.getGreen(), c.getBlue())
            }
            // Palette expects first color unused (index 0), histogram uses index i+1
            Palette p = new Palette()
            p.colors = ([colors[0]] + colors) as Color[]
            h.palette = p
        }

        return h
    }

    @CompileStatic
    private List<Double> calculateBreaks(Iterable<Iterable<Double>> datas) {
        
        Stats dataStats = new Stats()
        for(Iterable<Double> data in datas) {
            for(d in data) {
                if(d >= rangeMin && d <= rangeMax)
                    dataStats.addValue(d)
            }
        }

        double range = dataStats.max - dataStats.min
        double delta = (dataStats.max - dataStats.min + Double.MIN_VALUE) / (binCount-1);
        double halfDelta = delta/2d;

        List<Double> breaks = []
        double pos = dataStats.min - halfDelta
        for(int i=0; i<=binCount; ++i) {
            breaks.add(pos)
            pos += delta
        }
        return breaks
    }
}

class Plot {
    
    String title = 'Plot'
    
    String xLabel = null
    
    String yLabel = null
    
    List<PlotItem> items = []

    List xBound = null
    
    List yBound = null

    String legendLocation = null // String 
    
    List<Text> texts = []
    
    List<ConstantLine> constantLines = []
    
    Palette palette = new DefaultPalette()

    /**
     * Default style for tooltips displayed on this plot. Individual series may
     * override it via {@link XYItem#toolTipStyle}.
     */
    ToolTipStyle toolTipStyle = new ToolTipStyle()

    /**
     * For compatibility with BeakerX
     */
    Integer initWidth = null
    Integer initHeight = null
    
    Double legendDistance = null
    
    Plot leftShift(PlotItem item) {
        this.items << item
        return this
    }
    
    Plot leftShift(com.twosigma.beakerx.chart.xychart.plotitem.Area item) {
        return addBeakerXItem(item, new Area())
    }
     
    Plot leftShift(com.twosigma.beakerx.chart.xychart.plotitem.Line item) {
        return addBeakerXItem(item, new Line())
    }
     
    Plot leftShift(com.twosigma.beakerx.chart.xychart.plotitem.Points item) {
        return addBeakerXItem(item, new Points())
    }
     
    Plot leftShift(com.twosigma.beakerx.chart.xychart.plotitem.Bars item) {
        return addBeakerXItem(item, new Bars())
    }
     
    /**
     * Add the gngs equivalent of a BeakerX graphics item to this plot, carrying
     * across every property the two have in common.
     * 
     * @param source    the BeakerX item
     * @param item      the gngs item standing in for it, which is added to the plot
     */
    private Plot addBeakerXItem(XYGraphics source, XYItem item) {
        copyBeakerXProperties(source, item, this.items.size())
        this.items << item
        return this
    }
    
    /**
     * Copy the properties of a BeakerX chart object onto the gngs object that
     * stands in for it, matching them by name.
     * <p>
     * Copying reflectively rather than naming the fields keeps the two sides
     * from drifting apart as either gains attributes, which is how BeakerX line
     * styles came to be silently dropped. Tooltips are the one thing that has
     * to be carried over by hand, because the two sides do not agree on the name.
     * 
     * @param source    object to read properties from
     * @param item      object to copy them onto, ignoring any it does not have
     * @param index     position of the item, used to pick a palette colour
     */
    private void copyBeakerXProperties(Object source, Object item, int index) {
        
        source.properties.each { k, v ->
            
            if(k == "color" && v instanceof com.twosigma.beakerx.chart.Color) {
                v = convertColor(v, index)
            }

            if(item.hasProperty(k)) {
                try {
                    item[k] = v
                }
                catch(ReadOnlyPropertyException exReadOnly) {
                    // eg: the class property; nothing to be done and nothing wanted
                }
            }
        }

        // Not matched by name above: BeakerX exposes tooltips as the read only
        // property toolTips, while ours is called toolTip
        if(source instanceof XYGraphics && item instanceof XYItem) {
            Object tips = beakerXToolTips((XYGraphics)source)
            if(tips != null)
                ((XYItem)item).toolTip = tips
        }
    }
     
    Plot leftShift(com.twosigma.beakerx.chart.xychart.plotitem.Text item) {
        Text gngsItem = new Text(x:item.x, y:item.y, text: item.text)
        if(item.color)
            gngsItem.color = new Color(item.color.RGB)
        this.texts << gngsItem
        return this
    }

    Plot leftShift(Text item) {
        this.texts << item
        return this
    }
     
    Plot leftShift(com.twosigma.beakerx.chart.xychart.plotitem.ConstantLine item) {
        ConstantLine gngsItem = new ConstantLine(x:item.x, y:item.y, style: item.style)
        if(item.color)
            gngsItem.color = new Color(item.color.RGB)
        this.constantLines << gngsItem
        return this
    }

     Plot leftShift(ConstantLine item) {
        this.constantLines << item
        return this
    }

    /**
     * Display tooltips on every series of this plot that has them.
     * <p>
     * Example:
     * <pre>
     * def p = new Plot(title: 'An example of showing a tooltip') &lt;&lt;
     *     new Points(x: [1,2,3,4], y: [5,6,7,8], toolTip: [1,2,3,4].collect { "X value is: $it" })
     *
     * p.showTooltips([0,2])
     * p.save('plot.png')
     * </pre>
     *
     * @param style optional named arguments overriding {@link ToolTipStyle} properties
     * @param which which tooltips to show, see {@link XYItem#showTooltips(Map,Object)}
     */
    Plot showTooltips(Map style = null, Object which) {
        List<XYItem> annotatable = this.items.grep { it instanceof XYItem && it.toolTip != null } as List<XYItem>

        if(annotatable.isEmpty())
            throw new IllegalStateException('No series in this plot has any toolTip set')

        for(XYItem item in annotatable) {
            item.showTooltips(style?:[:], which)
        }
        return this
    }

    /**
     * Display tooltips on the series having the given display name
     *
     * @param style       optional named arguments overriding {@link ToolTipStyle} properties
     * @param displayName display name of the series to annotate
     * @param which       which tooltips to show, see {@link XYItem#showTooltips(Map,Object)}
     */
    Plot showTooltips(Map style = null, String displayName, Object which) {
        // Note: an explicit cast, not "as XYItem", because XYItem overrides
        // asType and returns null for anything that is not a DataTable
        XYItem item = (XYItem)this.items.find { it instanceof XYItem && it.displayName == displayName }

        if(item == null)
            throw new IllegalArgumentException(
                "No series found with display name '$displayName'. Available: " +
                (this.items*.displayName.grep { it } .join(', ') ?: '<none set>'))

        item.showTooltips(style?:[:], which)
        return this
    }

    /**
     * Display tooltips on the series at the given index, counting only the
     * items that are plotted as x/y data
     *
     * @param style  optional named arguments overriding {@link ToolTipStyle} properties
     * @param series index of the series to annotate
     * @param which  which tooltips to show, see {@link XYItem#showTooltips(Map,Object)}
     */
    Plot showTooltips(Map style = null, int series, Object which) {
        List<XYItem> xys = this.items.grep { it instanceof XYItem } as List<XYItem>

        if(series < 0 || series >= xys.size())
            throw new IllegalArgumentException("Series $series does not exist: this plot has ${xys.size()} series")

        xys[series].showTooltips(style?:[:], which)
        return this
    }

    /**
     * Display every tooltip set on every series of this plot
     *
     * @param style optional named arguments overriding {@link ToolTipStyle} properties
     */
    Plot showAllTooltips(Map style = null) {
        for(XYItem item in this.items.grep { it instanceof XYItem && it.toolTip != null }) {
            item.showAllTooltips(style?:[:])
        }
        return this
    }

    static Object saveAs(def plot, String fileName) {
        if(plot instanceof com.twosigma.beakerx.chart.xychart.Plot) {
            from(plot).save(fileName)
        }
        else
        if(plot instanceof com.twosigma.beakerx.chart.histogram.Histogram) {
            Histogram.from(plot).save(fileName)
        }
        else
        if(plot instanceof Plot) {
            plot.save(fileName)
        }
        else
        if(plot instanceof Histogram) {
            plot.save(fileName)
        }
        else {
            throw new IllegalArgumentException('Please provide a gngs or beakerx Plot or Histogram object - you provided: ' + plot?.class?.name)
        }
        return plot
    }
    
    static from(gngs.plot.Plot other) {
        return other
    }

    static from(gngs.plot.Histogram other) {
        return other
    }

    static from(com.twosigma.beakerx.chart.histogram.Histogram bxHist) {
        return Histogram.from(bxHist)
    }

    static Plot from(com.twosigma.beakerx.chart.xychart.Plot bxPlot) {
        Plot p = new Plot(
            title:bxPlot.title,
            xLabel: bxPlot.xLabel,
            yLabel: bxPlot.yLabel,
            initWidth: bxPlot.initWidth,
            initHeight: bxPlot.initHeight
        )
        
        if(!bxPlot.xAutoRange)
            p.xBound = [bxPlot.xLowerBound, bxPlot.xUpperBound]

        if(!bxPlot.yAutoRange)
            p.yBound = [bxPlot.getYLowerBound(), bxPlot.getYUpperBound()]
        
        def setProps = { g, item, i ->
            p.copyBeakerXProperties(g, item, (int)i)
        }

        int i = 0
        bxPlot.graphics.each { XYGraphics g ->
            
            def item = null
            
            if(g instanceof com.twosigma.beakerx.chart.xychart.plotitem.Points) {
                item = new Points()
            }
            else
            if(g instanceof com.twosigma.beakerx.chart.xychart.plotitem.Line) {
                item = new Lines()
            }
            else
            if(g instanceof com.twosigma.beakerx.chart.xychart.plotitem.Area) {
                item = new Area()
            }
            if(g instanceof com.twosigma.beakerx.chart.xychart.plotitem.Bars) {
                item = new Bars()
            }
            if(g instanceof com.twosigma.beakerx.chart.xychart.plotitem.Text) {
                item = new Text()
            }
     
            if(!item)
                return

            setProps(g, item, i)

            p << item
            ++i
        }
        
        bxPlot.constantLines.each { cl ->
            ConstantLine item = new ConstantLine()
            setProps(cl, item, i)
            p << item
            ++i
        }
        
        return p
    }
    
    /**
     * Extract the tooltips from a BeakerX graphics item.
     * <p>
     * BeakerX accepts tooltips either as a literal list or as a closure which
     * builds the text for each point, but it resolves the closure eagerly at the
     * point it is assigned, so there is always a resolved list to read here.
     *
     * @return the tooltips, or null if the item has none
     */
    private static Object beakerXToolTips(XYGraphics g) {
        List<String> tips = g.getToolTips()
        return tips ?: null
    }

    BufferedImage getImage() {
        getImage(initWidth?:800,initHeight?:600)
    }

    BufferedImage getImage(int width, int height) {

        int eastLegendWidth = estimateLegendWidth(width)

        int rasterFormat = BufferedImage.TYPE_INT_RGB;
        BufferedImage image = new BufferedImage(
                (int)Math.ceil(width + eastLegendWidth), (int)Math.ceil(height), rasterFormat);

        DrawingContext context = PlotUtils.createDrawingContext(image)

        XYPlot xyPlot = toXYPlot(width,height)
        
        Rectangle2D boundsOld = xyPlot.getBounds();
        xyPlot.setBounds(0, 0, width, height);
        xyPlot.draw(context)
        
        return image
    }
    
    XYPlot toXYPlot(int imageWidth, int imageHeight) {
       
        List<XYItem> xys = items.grep { it instanceof XYItem }
        
        int uberMinX = xys*.maxX.min()
        int uberMaxX = PlotUtils.roundUpToOOM(xys*.maxX.max())
        int uberMinY = xys*.minY.min()
        int uberMaxY = PlotUtils.roundUpToOOM(xys*.maxY.max())
        
        int i = 1
        
        List<DataTable> clDatas = constantLines.collect { ConstantLine cl ->
            DataTable dt  = cl.toTable(uberMinX, uberMaxX, uberMinY, uberMaxY) 
            if(cl.displayName) {
                dt.setName(cl.displayName)
            }
            return dt
         }
        
        List<DataTable> datas = xys.collect { XYItem item ->
            DataTable dt = item.toTable()
            dt.setName(item.displayName ?: ('Series ' + i))
//            dt.setName(item.displayName)
            ++i
            return dt
        } + clDatas
        
        DataTable [] dtArray = datas as DataTable[]
        
        XYPlot xyPlot = 
            xys.any { it instanceof Bars } ? 
                new BarPlot(dtArray)
            :
                new XYPlot(dtArray)
                
        XYPlot legendPlot = new XYPlot(dtArray)

        Insets2D.Double insets = new Insets2D.Double(40.0, 80.0, 80.0, 80.0)
        xyPlot.setInsets(insets);
        
        i = 0
        
        [xys,datas].transpose().each { xy, dt ->

            def color = convertColor(xy.color, i)

            if(xy instanceof Lines || xy instanceof Line) {
                LineRenderer lines = new SmoothLineRenderer2D();

                Stroke stroke = null
                if(xy.width != null)
                    stroke = new BasicStroke((float)xy.width)

                if(xy.style == "DOT") {
                    stroke = new BasicStroke(
                        xy.width ? xy.width.toFloat() :1.0f,           // line width
                        BasicStroke.CAP_ROUND,    // round caps for dot effect
                        BasicStroke.JOIN_MITER,   // round joins
                        4.0f,                     // miter limit
                        [1f, 6f] as float[],      // pattern: 1px dash, 4px gap
                        0.0f                      // phase offset
                    );
                }
                else
                if(xy.style == "DASH") {
                    stroke = new BasicStroke(
                        xy.width ? xy.width.toFloat() :1.0f,           // line width
                        BasicStroke.CAP_ROUND,    // round caps for dot effect
                        BasicStroke.JOIN_MITER,   // round joins
                        4.0f,                     // miter limit
                        [8f, 8f] as float[],      // pattern: 1px dash, 4px gap
                        0.0f                      // phase offset
                    );
                }
                    

                if(stroke)
                    lines.setStroke(stroke)

                lines.setColor(color)
                xyPlot.setLineRenderers(dt, lines)
                xyPlot.setPointRenderers(dt, null)
            }
            else
            if(xy instanceof Area) {

                AreaRenderer area = new DefaultAreaRenderer2D();
                area.color = GraphicsUtils.deriveWithAlpha(color, 115)
                area.setGapRounded(true)
                
                PointRenderer point = new DefaultPointRenderer2D();
                point.setColor(area.color);
                xyPlot.setPointRenderers(dt, point);
                    
                
                LineRenderer line = new DefaultLineRenderer2D();
                line.setColor(xy.color);
                line.setGap(0);
                xyPlot.setLineRenderers(dt, line);
                xyPlot.setAreaRenderers(dt, area)
            }
            else
            if(xy instanceof Bars) {
                List<BarRenderer> bars = xyPlot.getPointRenderers(dtArray[i])
                bars*.setColor(color)
                if(xy.width != null)
                    xyPlot.setBarWidth(xy.width)
                    
                legendPlot.getPointRenderers(dtArray[i])*.setColor(color)
                    
                if(xy.labels != null) {
                    bars*.setValueColumn(2)
                    bars*.setValueVisible(true)
                }
            }
            else {
                PointRenderer pointRenderer = new DefaultPointRenderer2D()
                pointRenderer.setColor(palette.colors[ i % palette.colors.size()])
                xyPlot.setPointRenderers(dt, pointRenderer)
                pointRenderer.setColor(color)

            }
//            xyPlot.setMapping(dt, xys[i].name, '')
            ++i
        }
        
//        for(dt in datas) {
//            LineRenderer lines = new SmoothLineRenderer2D();
//            lines.setColor(palette.colors[ i % palette.colors.size()])
//            xyPlot.setLineRenderers(dt, lines)
////            xyPlot.setMapping(dt, xys[i].name, '')
//            ++i
//        }
//        
        xyPlot.setBackground(Color.white)
        xyPlot.getTitle().setText(title)        
        
        def xAxis = xyPlot.getAxis(XYPlot.AXIS_X)
        xAxis.with {
            if(xBound) {
                min = xBound[0]
                max = xBound[1]
            }
            else {
                min = Math.min(0, xys*.minX.min())
                max = PlotUtils.roundUpToOOM(xys*.maxX.max())
            }
        }
        
        def yAxis = xyPlot.getAxis(XYPlot.AXIS_Y)
        yAxis.with {
            
            if(yBound) {
                min = yBound[0]
                max = yBound[1]
            }
            else {
                min = Math.min(0, xys*.minY.min())
                max = PlotUtils.roundUpToOOM(xys*.maxY.max())
            }
        }
        
        xyPlot.getAxisRenderer(XYPlot.AXIS_X).with { 
            if(yBound) {
                intersection = (yBound[0] as double)
            }

            if(xLabel)
                label.text = xLabel
        }
        
        xyPlot.getAxisRenderer(XYPlot.AXIS_Y).with { 
            if(xBound) {
                intersection = (xBound[0] as double)
            }
            if(yLabel)
                label.text = yLabel
        }
        
        if(xys.any { it.displayName }) {
            
            if(xyPlot instanceof BarPlot) {
                // noop
                // default legend doesn't work here because it displays
                // a legend entry for each bar rather than each series
                // which is the default for bar plots.
                // Need to replace with a SeriesLegend, but it seems like
                
                if(xys.size() > 1) {
                    SeriesLegend legend = new XYLegend(legendPlot)
                    xyPlot.setLegend(legend)
                    xyPlot.setLegendVisible(true)
                    if(this.legendLocation) {
                        xyPlot.setLegendLocation(Location[this.legendLocation.toUpperCase()])
                    }
                }
            }
            else {
                xyPlot.setLegendVisible(true)
                
                if(this.legendLocation) {
                    xyPlot.setLegendLocation(Location[this.legendLocation.toUpperCase()])
                }
            }
            if(this.legendDistance != null) {
                xyPlot.setLegendDistance(this.legendDistance)
            }
        }
            
        double width = (double)imageWidth
        double height = (double)imageHeight
        addTextsToXY(xyPlot, xAxis, yAxis, width, height)
        
        clDatas.each {
            addConstantLinesToXYPlot(it, xyPlot, xAxis, yAxis, width, height)
        }

        addToolTipsToXYPlot(xyPlot, xys)

        return xyPlot
    }

    /**
     * Add an annotation layer for any tooltips that have been selected for display.
     * <p>
     * The layer resolves its own positions when it is drawn, because the plot
     * has not been laid out at the point where it is added here.
     */
    void addToolTipsToXYPlot(XYPlot xyPlot, List<XYItem> xys) {

        List<ToolTipAnnotation> annotations = []

        xys.eachWithIndex { XYItem item, int series ->

            if(!item.shownToolTips)
                return

            List<String> tips = item.resolveToolTips()
            ToolTipStyle style = item.toolTipStyle ?: this.toolTipStyle

            item.shownToolTips.each { Integer index, ToolTipPlacement placement ->

                if(index < 0 || index >= tips.size())
                    throw new IllegalArgumentException(
                        "Tooltip index $index is out of range for series " +
                        "'${item.displayName?:series}' which has ${tips.size()} tooltips")

                String text = tips[index]
                if(!text)
                    return

                annotations << new ToolTipAnnotation(
                    series: series,
                    index: index,
                    text: text,
                    placement: placement,
                    style: style
                )
            }
        }

        if(annotations.isEmpty())
            return

        xyPlot.add(new ToolTipLayer(plot: xyPlot, items: xys, annotations: annotations))
    }
    
    void addConstantLinesToXYPlot(DataTable dt, XYPlot xyPlot, Axis xAxis, Axis yAxis, double width, double height) {

        for(ConstantLine cl in constantLines) {
            // Add a constant line as if it was a two point line
            LineRenderer lines = new SmoothLineRenderer2D();
//            LineRenderer lines = new DiscreteLineRenderer2D();
            
            BasicStroke stroke = null
            // Unfortunately setting the stroke with a pattern causes an out of memory error - really not sure why, but
            // it could interact badly with the fact we only have two points in our plot maybe?
//            if(cl.style == "DOT") {
////                float[] dashPattern = [1f, 1f] as float[]; // small dots with spacing
////                stroke = new BasicStroke(1f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_MITER, 2f, dashPattern, 0f);
//                
//                stroke = new BasicStroke(
//                    1.0f,                      // line width
//                    BasicStroke.CAP_ROUND,    // round caps for dot effect
//                    BasicStroke.JOIN_MITER,   // round joins
//                    4.0f,                     // miter limit
//                    new float[]{1f, 4f},      // pattern: 1px dash, 4px gap
//                    0.0f                      // phase offset
//                );
//            }
//            else
            if(cl.width != null)
                stroke = new BasicStroke((float)cl.width)
                
            if(stroke) {
                lines.setStroke(stroke)
            }

            def color = convertColor(cl.color, 0)
            lines.setColor(color)
                
            xyPlot.setLineRenderers(dt, lines)
        }
    }
    
    void addTextsToXY(XYPlot xyPlot, Axis xAxis, Axis yAxis, double width, double height) {
        Insets2D.Double insets = xyPlot.getInsets()
        for(Text text in texts) {

            Label label = new Label(text.text)
           
            int renderXOffset = xyPlot.getAxisRenderer(XYPlot.AXIS_X).worldToView(xAxis, 0.0d, false)
            int labelRenderX = insets.left + width * (text.x - xAxis.min) / (xAxis.max - xAxis.min) - renderXOffset
            int labelRenderY = height - (insets.top + width * (text.y - yAxis.min) / (yAxis.max - yAxis.min))
            
            label.setPosition(labelRenderX,labelRenderY)
            label.alignmentX=0.0d
            label.background = Color.orange
            
            if(text.color) {
                label.color = text.color
            }

            xyPlot.add(label)
        }
    }
    
    /**
     * Convert to Gral compatible object, using palette for item index if not
     * explicitly specified
     * 
     * @param color color object to convert, may be null
     * @param i index to select from palette if color is null
     * @return  Gral compatible color object
     */
    Color convertColor(def color, int i) {
        
        if(color)
            // This slightly awkward way works regardless of the type of color object passed in
            // ie: AWT, BeakerX, etc
            return new Color(color.red, color.green, color.blue)
        else
            return palette.colors[ i % palette.colors.size()]
    }
    
    @CompileStatic
    void save(Map options=null, final String fileName) {
        
        assert fileName.endsWith('.png')

        options = options?:[:]

        int width = (int)(options.width?:initWidth?:1024)
        int height = (int)(options.height?:initHeight?:800)
        
        XYPlot xyPlot = toXYPlot(width, height)
        
        int eastLegendWidth = (int)(options.marginRight ?: this.estimateLegendWidth(width))

        new File(fileName).withOutputStream { w ->
            DrawableWriter wr = DrawableWriterFactory.getInstance().get("image/png");
            PlotUtils.write(xyPlot, w, 0,0, width, height, eastLegendWidth);
        } 
    }
    
    @CompileStatic
    int estimateLegendWidth(int width) {
        int eastLegendWidth = 0
        if(this.legendLocation in ["east","north_east","south_east"]) {
            int maxDisplayNameLength = this.items*.displayName.collect { it?.size()?:0 }.max()
            eastLegendWidth = (int)(maxDisplayNameLength*10 * (width/1024)) // hack / guess, use about 15% of width
        }
        return eastLegendWidth
    }
}
    
class PlotUtils {
    
    @CompileStatic
    static double roundUpToOOM(double x) {
        if(x == 0d)
            return 1d
        int oom = (int)Math.floor(Math.log10(x)) 
        double interval = Math.pow(10,oom)
        double rounded = Math.floor(x/interval) * interval  + interval
        return rounded
    }
    
    /*
     * Stores the specified {@code Drawable} instance.
     * @param d {@code Drawable} to be written.
     * @param destination Stream to write to
     * @param x Horizontal position.
     * @param y Vertical position.
     * @param width Width of the image.
     * @param height Height of the image.
     * @throws IOException if writing to stream fails
     */
    static public void write(Drawable d, OutputStream destination,
            double x, double y, double width, double height, int legendWidth, Closure then = null)
            throws IOException {

        int rasterFormat = BufferedImage.TYPE_INT_RGB;
        BufferedImage image = new BufferedImage(
                (int)Math.ceil(width + legendWidth), (int)Math.ceil(height), rasterFormat);

        DrawingContext context = createDrawingContext(image)

        Iterator<ImageWriter> writers =
            ImageIO.getImageWritersByMIMEType('image/png');
        if (writers.hasNext()) {
            ImageWriter writer = writers.next();
            ImageOutputStream ios =
                ImageIO.createImageOutputStream(destination);
            writer.setOutput(ios);
            Rectangle2D boundsOld = d.getBounds();
            d.setBounds(x, y, width, height);
            try {
                d.draw(context);
                if(then) {
                    then(context)
                }
                writer.write(image);
            } finally {
                d.setBounds(boundsOld);
                ios.close();
            }
        }
    }
    
    static DrawingContext createDrawingContext(BufferedImage image) {
        Graphics2D imageGraphics = image.createGraphics();
        imageGraphics.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        imageGraphics.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
        imageGraphics.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
        imageGraphics.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_BICUBIC);
        imageGraphics.setRenderingHint(RenderingHints.KEY_FRACTIONALMETRICS, RenderingHints.VALUE_FRACTIONALMETRICS_ON);
        
        imageGraphics.background = Color.white
        imageGraphics.fillRect(0, 0, image.width, image.height)
        
        DrawingContext context = new DrawingContext(imageGraphics);        
        
        return context
    }
}
