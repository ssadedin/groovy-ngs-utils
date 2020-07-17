package gngs.plot.bx

import com.twosigma.beakerx.chart.Color
import com.twosigma.beakerx.chart.xychart.plotitem.XYGraphics
import groovy.transform.CompileStatic

@CompileStatic
class DensityPoint {
    public double x
    public double y
}

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
    
    smile.stat.distribution.KernelDensity kd
    
    double step
    
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
        
        def me = this
        
        attributes.each { String k, Object v -> 
            if(!(k in DENSITY_ATTRIBUTES))  {
                if(k == 'color' && (v instanceof Color)) {
                    me.setColor(v)
                }
                else {
                    String method = 'set' + k.capitalize()
                    me.invokeMethod(method,[v])
                }
            }
        }
        
        com.twosigma.beakerx.chart.xychart.plotitem.Area area
        
        if(!('data' in attributes)) 
            throw new IllegalArgumentException('Please provide a data attribute with values to display')
            
        double [] values = attributes.data as double[]
        double min = ((List)attributes.data).min()
        double max = ((List)attributes.data).max()
        
        initXYValues(min, max, values)
    }    
    
    @CompileStatic
    void initXYValues(double min, double max, double[] values) {
        
        kd = new smile.stat.distribution.KernelDensity(values, bw)
        step = (max - min) / segments
        
        min -= step*2
        max += step*2
        
        final double plotMin = min - step;
        final double plotMax = max + step;
        
        LinkedList<DensityPoint> plotValues = new LinkedList()

        double maxDensity = computeCoreValues(plotMin, plotMax, plotValues)

        expandLeading(plotMin, maxDensity, plotValues)

        expandTrailing(plotMax, maxDensity, plotValues)

        this.x = plotValues*.x
        this.y = plotValues*.y        
    }
    
    @CompileStatic
    double computeCoreValues(double plotMin, double plotMax, List plotValues) {
        double maxDensity = 0.0d
        for(double value = plotMin; value < plotMax; value += step) {
            double density = kd.p(value)
            plotValues <<  new DensityPoint(y:density, x:value)
            if(density > maxDensity)
                maxDensity = density
        }
        return maxDensity
    }
    
    @CompileStatic
    void expandTrailing(final double plotMax, final double maxDensity, final List plotValues) {
        double tailing = plotMax + step
        while(true) {
            double density = kd.p(tailing)
            if(density < 0.05 * maxDensity) {
                break
            }
            plotValues.add(new DensityPoint(y:density, x:tailing))
            tailing += step
        }
    }
    
    @CompileStatic
    void expandLeading(final double plotMin, final double maxDensity, final List plotValues) {

        // Problem: density plots usually spill pass the edge of the sampled data points
        // therefore we keep steping backwards until the tail is below 5% of the max
        // which ensures the distribution does not get too clipped at the end
        double leading = plotMin - step
        while(true) {
            double density = kd.p(leading)
            if(density < 0.05 * maxDensity) {
                break
            }
            plotValues.add(0, new DensityPoint(y:density, x:leading))
            leading -= step
        }
    }
 }
