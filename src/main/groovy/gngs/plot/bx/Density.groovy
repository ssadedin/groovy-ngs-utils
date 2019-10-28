package gngs.plot.bx

import com.twosigma.beakerx.chart.xychart.plotitem.XYGraphics

/**
 * Provides support for plotting of empirical PDFs in BeakerX
 * <p>
 * Plotting is implemented as specific subclasses for each of the different
 * high level plotting elements (Line,Area,Points).
 * <p>
 * These operate the same way that regular usage of these elements works, however
 * instead of providing x and y attributes, the data values are provided for which
 * the probability density should be estimated, and the x and y values are derived 
 * from those.
 * <p>
 * Example:
 * <pre>
 * def x = (1..1000).collect { r.nextGaussian()*20 }
 * new Plot(title:'An example plot') << new Density.Line(data:x, color: Color.red)
 * </pre>
 * @author Simon Sadedin
 */
trait Density {
    
    final static List DENSITY_ATTRIBUTES = ['data','segments','bw']
    
    double bw = 2.0d
    
    int segments = 100
    
    static class Line extends  com.twosigma.beakerx.chart.xychart.plotitem.Line implements Density {
        Line(Map attributes) { init(attributes) }
    }
    
    static class Area extends  com.twosigma.beakerx.chart.xychart.plotitem.Area implements Density {
        Area(Map attributes) { init(attributes) }
    }
    
    static class Points extends  com.twosigma.beakerx.chart.xychart.plotitem.Points implements Density {
        Points(Map attributes) { 
            init(attributes) ;
        }
    } 
    
    void init(Map attributes) {
        
        this.bw = attributes.get('bw',bw).toDouble()
        this.segments = attributes.get('segments',segments).toInteger()
        
        attributes.each { k,v -> if(!(k in DENSITY_ATTRIBUTES)) this[k] = v; }
        
        if(!('data' in attributes)) 
            throw new IllegalArgumentException('Please provide a data attribute with values to disply')
            
        double [] values = attributes.data as double[]
        def kd = new smile.stat.distribution.KernelDensity(values, bw)
        double min = attributes.data.min()
        double max = attributes.data.max()
        double step = (max - min) / segments
        
        min -= step
        max += step
        
        List yValues = new ArrayList(segments)
        List xValues = new ArrayList(segments)
        for(double value = min; value < max; value += step) {
            yValues << kd.p(value)
            xValues << value
        }
        
        this.x = xValues
        this.y = yValues
    }    
 }
