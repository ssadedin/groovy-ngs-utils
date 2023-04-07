package gngs.plot.bx

import org.apache.commons.math3.stat.StatUtils

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
    
    Double bw = null
    
    double cutoff = 0.01
    
    boolean cumulative = false
    
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
        
        this.bw = attributes.get('bw',bw)?.toDouble()
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
            
        Iterable valuesIterable
        double [] values 
        double min
        double max
        if(attributes.data instanceof double[]) {
            values = attributes.data
            max = StatUtils.max(values)
            min = StatUtils.min(values)
        }
        else {
            values = attributes.data as double[]
            min = ((Iterable)attributes.data).min()
            max = ((Iterable)attributes.data).max()
        }
        
        initXYValues(min, max, values)
    }    
    
    @CompileStatic
    void initXYValues(double min, double max, double[] values) {
        
        kd = (bw != null ? 
            new smile.stat.distribution.KernelDensity(values, bw) : 
            new smile.stat.distribution.KernelDensity(values)) 
             
        step = (max - min) / segments
        
        // This is a bit arbitrary, but if all the values are the same, we will just show
        // a bin of 1.0 either side of whatever the value is
        if(step == 0.0d) {
            step = 1.0d
        }
        
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
            double density = cumulative ? kd.cdf(value) : kd.p(value)
            plotValues <<  new DensityPoint(y:density, x:value)
            if(density > maxDensity)
                maxDensity = density
        }
        return maxDensity
    }
    
    final int maxZeroGradientSteps = 10
    
    @CompileStatic
    void expandTrailing(final double plotMax, final double maxDensity, final List plotValues) {
        double tailing = plotMax + step

        // Need to handle problem when gradient is zero: in that case
        // we will never converge to a value that satisfies the criteria
        int zeroGradientSteps = 0
        double lastDensity = Double.POSITIVE_INFINITY

        while(true) {
            double density = cumulative ? kd.cdf(tailing) : kd.p(tailing)
            if(density-lastDensity == 0.0d) {
                ++zeroGradientSteps
                if(zeroGradientSteps > maxZeroGradientSteps)
                    break
            }
            else {
                zeroGradientSteps = 0
            }
            lastDensity = density

            if(density < cutoff * maxDensity) {
                break
            }
            plotValues.add(new DensityPoint(y:density, x:tailing))
            tailing += step
        }
    }
    
    @CompileStatic
    void expandLeading(final double plotMin, final double maxDensity, final List plotValues) {

        // Problem: density plots usually spill pass the edge of the sampled data points
        // therefore we keep stepping backwards until the tail is below 5% of the max
        // which ensures the distribution does not get too clipped at the end
        double leading = plotMin - step
        
        // Need to handle problem when gradient is zero: in that case
        // we will never converge to a value that satisfies the criteria
        int zeroGradientSteps = 0
        double lastDensity = Double.POSITIVE_INFINITY
        while(true) {
            double density = cumulative ? kd.cdf(leading) : kd.p(leading)
            if(density-lastDensity == 0.0d) {
                ++zeroGradientSteps
                if(zeroGradientSteps > maxZeroGradientSteps)
                    break
            }
            else {
                zeroGradientSteps = 0
            }
            lastDensity = density
            
            if(density < cutoff * maxDensity) {
                break
            }
            plotValues.add(0, new DensityPoint(y:density, x:leading))
            leading -= step
        }
    }
 }
