package gngs.plot

import java.awt.Color
import java.awt.Graphics2D
import java.awt.RenderingHints
import java.awt.geom.Rectangle2D
import java.awt.image.BufferedImage

import javax.imageio.ImageIO
import javax.imageio.ImageWriter
import javax.imageio.stream.ImageOutputStream
import de.erichseifert.gral.data.DataSource
import de.erichseifert.gral.data.DataTable
import de.erichseifert.gral.data.EnumeratedData
import de.erichseifert.gral.data.statistics.Histogram2D
import de.erichseifert.gral.graphics.Drawable
import de.erichseifert.gral.graphics.DrawingContext
import de.erichseifert.gral.graphics.DrawingContext.Quality
import de.erichseifert.gral.graphics.DrawingContext.Target
import de.erichseifert.gral.graphics.Insets2D
import de.erichseifert.gral.graphics.Label
import de.erichseifert.gral.graphics.Orientation
import de.erichseifert.gral.io.plots.DrawableWriter
import de.erichseifert.gral.io.plots.DrawableWriterFactory
import de.erichseifert.gral.plots.BarPlot
import de.erichseifert.gral.plots.XYPlot
import de.erichseifert.gral.plots.axes.AxisRenderer
import de.erichseifert.gral.plots.lines.LineRenderer
import de.erichseifert.gral.plots.lines.SmoothLineRenderer2D
import de.erichseifert.gral.plots.points.PointRenderer
import de.erichseifert.gral.util.GraphicsUtils
import graxxia.Stats
import groovy.transform.CompileStatic


class Palette {
    Color [] colors 
    
}

class DefaultPalette extends Palette {
    DefaultPalette() {
       colors = [Color.black, Color.red, Color.blue, Color.green, Color.orange, new Color(100,0,100)] 
    }
}

class PlotItem {
    String name = null
}

class XYItem extends PlotItem {
    Iterable<Number> x
    Iterable<Number> y
    
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
        DataTable dt = new DataTable(2, Double)
        [x,y].transpose().each { Number xVal, Number yVal ->
            if(xVal > maxX)
                maxX = xVal
            if(yVal > maxY)
                maxY = yVal
            dt.add(xVal.toDouble(), yVal.toDouble())
        }
        return dt
    }
}

class Lines extends XYItem {
    
}

class Bars extends XYItem {
    
}

class Histogram {
    
    String title
    
    Iterable<Double> data
    
    int binCount
    
    Palette palette = new DefaultPalette()
    
    String xLabel
    
    String yLabel
    
    void save(final String fileName) {
        
        DataTable dt = new DataTable(1, Double)
        for(double d in data) {
            dt.add(d)
        }
        
        assert dt.columnCount == 1
        
        List breaks = calculateBreaks()
        
        double width = breaks[1] - breaks[0]
        
        // Create histogram from data
        Histogram2D hist = new Histogram2D(dt, Orientation.VERTICAL, [breaks as Double[]] as Double[][]);
        
        // Create a second dimension (x axis) for plotting
        DataSource histogram2d = new EnumeratedData(hist, breaks.min()+width/2, width);

        // Create new bar plot
        BarPlot plot = new BarPlot(histogram2d);
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

        plot.getAxis(XYPlot.AXIS_Y).with {
            max = PlotUtils.roundUpToOOM(hist.getColumn(0).max())
        }
        
        PointRenderer barRenderer = plot.getPointRenderers(histogram2d).get(0);
        barRenderer.setColor(GraphicsUtils.deriveWithAlpha(palette.colors[1], 128));
//        barRenderer.setValueVisible(true);
        
        Insets2D.Double insets = new Insets2D.Double(40.0, 80.0, 80.0, 80.0)
        plot.setInsets(insets);
        
        BufferedImage bImage = new BufferedImage(1024, 800, BufferedImage.TYPE_INT_ARGB);
        Graphics2D g2d = (Graphics2D) bImage.getGraphics();
        DrawingContext context = new DrawingContext(g2d, Quality.QUALITY, Target.BITMAP);
        
        new File(fileName).withOutputStream { w ->
            DrawableWriter wr = DrawableWriterFactory.getInstance().get("image/png");
            PlotUtils.write(plot, w, 0,0, 1024, 800);
        }        
    }

    @CompileStatic
    private List<Double> calculateBreaks() {
        Stats dataStats = Stats.from(data)
        double range = dataStats.max - dataStats.min
        double delta = (dataStats.max - dataStats.min + Double.MIN_VALUE) / (binCount-1);
        double halfDelta = delta/2d;

        List breaks = []
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
    
    Palette palette = new DefaultPalette()
    
    Plot leftShift(PlotItem item) {
        this.items << item
        return this
    }
    
    void save(final String fileName) {
        
        assert fileName.endsWith('.png')
        
        double maxX = Double.MIN_VALUE
        double maxY = Double.MIN_VALUE
        
        List<XYItem> xys = items.grep { it instanceof XYItem }
        
        int i = 1
        List<DataTable> datas = xys.collect { XYItem item ->
            DataTable dt = item.toTable()
            dt.name = item.name ?: ('Series ' + i)
            dt.setName(item.name)
            ++i
            return dt
        }
        
        XYPlot xyPlot = 
            xys.any { it instanceof Bars } ? 
                new BarPlot(datas as DataTable[])
            :
                new XYPlot(datas as DataTable[])
                
        Insets2D.Double insets = new Insets2D.Double(40.0, 80.0, 80.0, 80.0)
        xyPlot.setInsets(insets);
        
        i = 0
        for(dt in datas) {
            LineRenderer lines = new SmoothLineRenderer2D();
            lines.setColor(palette.colors[ i % palette.colors.size()])
            xyPlot.setLineRenderers(dt, lines)
//            xyPlot.setMapping(dt, xys[i].name, '')
            ++i
        }
        
        xyPlot.setBackground(Color.white)
        xyPlot.getTitle().setText(title)        
        
        xyPlot.getAxis(XYPlot.AXIS_X).with {
            min = Math.min(0, xys*.minX.min())
            max = PlotUtils.roundUpToOOM(xys*.maxX.max())
        }
        
        xyPlot.getAxis(XYPlot.AXIS_Y).with {
            min = Math.min(0, xys*.minY.min())
            max = PlotUtils.roundUpToOOM(xys*.maxY.max())
        }
        
        xyPlot.getAxisRenderer(XYPlot.AXIS_X).with { 
            if(xLabel)
                label.text = xLabel
        }
        
        xyPlot.getAxisRenderer(XYPlot.AXIS_Y).with { 
            if(yLabel)
                label.text = yLabel
        }
        
        if(!(xyPlot instanceof BarPlot))
            xyPlot.setLegendVisible(true)
  
        BufferedImage bImage = new BufferedImage(1024, 800, BufferedImage.TYPE_INT_ARGB);
        Graphics2D g2d = (Graphics2D) bImage.getGraphics();
        DrawingContext context = new DrawingContext(g2d, Quality.QUALITY, Target.BITMAP);
        
        new File(fileName).withOutputStream { w ->
            DrawableWriter wr = DrawableWriterFactory.getInstance().get("image/png");
            PlotUtils.write(xyPlot, w, 0,0, 1024, 800);
        } 
    }
    
}
    
class PlotUtils {
    
    @CompileStatic
    static double roundUpToOOM(double x) {
        if(x == 0d)
            return 1d
        int oom = (int)Math.floor(Math.log10(x)) 
        int interval = (int)Math.pow(10,oom)
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
            double x, double y, double width, double height)
            throws IOException {
        int rasterFormat = BufferedImage.TYPE_INT_RGB;
        BufferedImage image = new BufferedImage(
                (int)Math.ceil(width), (int)Math.ceil(height), rasterFormat);
        Graphics2D imageGraphics = image.createGraphics();
        imageGraphics.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        imageGraphics.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
        imageGraphics.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
        imageGraphics.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_BICUBIC);
        imageGraphics.setRenderingHint(RenderingHints.KEY_FRACTIONALMETRICS, RenderingHints.VALUE_FRACTIONALMETRICS_ON);
        

        DrawingContext context =
            new DrawingContext(imageGraphics);

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
                writer.write(image);
            } finally {
                d.setBounds(boundsOld);
                ios.close();
            }
        }
    }
}
