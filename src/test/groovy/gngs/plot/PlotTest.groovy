package gngs.plot

import static org.junit.Assert.*

import gngs.plot.bx.Density
import java.awt.Color

import org.junit.Test

class PlotTest {

//    @Test
//    public void testLines() {
//        Plot p = new Plot(
//            title:'Simons great plot',
//            xLabel: "Random Stuff You Don't Care About",
//            yLabel: "Stuff You Do Care About",
//            ) << \
//            new Lines(x: [1,2,3,4,5,6,7], y: [1,3,6,8,6,5,2], name: 'Bananas') << \
//            new Lines(x: [1,2,3,4,5,6,7], y: [1,5,9,18,4,2,1], name: 'Oranges')
//        p.save('test.png')
//    }
//    
//   @Test
    public void testPoints() {
        Plot p = new Plot(
            title:'Simons great plot',
            xLabel: "Random Stuff You Don't Care About",
            yLabel: "Stuff You Do Care About",
            xBound: [-20,20],
            yBound: [-5, 30]
            ) << \
            new Points(x: [1,2,3,4,5,6,7], y: [1,3,6,8,6,5,2], displayName: 'Bananas') << \
            new Points(x: [1,2,3,4,5,6,7], y: [1,5,9,18,4,2,1], displayName: 'Oranges') << \
            new Text(x: 0, y:0, text: 'Hello there world', color: Color.green)

        p.save('testpoints.png')
    }
  
//    @Test // fails outside beakerx :-(
    public void testBeakerx() {
        def p = new com.twosigma.beakerx.chart.xychart.Plot(
            title:'Simons great plot',
            xLabel: "Random Stuff You Don't Care About",
            yLabel: "Stuff You Do Care About",
            xBound: [-20,20],
            yBound: [-5, 30]
            ) << \
            new com.twosigma.beakerx.chart.xychart.plotitem.Points(x: [1,2,3,4,5,6,7], y: [1,3,6,8,6,5,2], displayName: 'Bananas')

        Plot gngsPlot = Plot.from(p)

//        p.save('testpoints.png')
    }
    
    @Test
    public void testBars() {
        Plot p = new Plot(title:'Simons great plot') << \
            new Bars(x: [1,2,3,4,5,6,7], y: [1,3,6,8,6,5,2], displayName: 'Bananas') << \
            new Bars(x: [1,2,3,4,5,6,7], y: [1,5,9,18,4,2,1], displayName: 'Oranges')
        p.save('testbars.png')        
    }
    
    
    
//    @Test
    void 'test round up oom' () {
        
        def p = new Plot()
        
        assert PlotUtils.roundUpToOOM(7d) == 10d
        
        assert PlotUtils.roundUpToOOM(17d) == 20d
        
    }
    
    
//    @Test
    void testBarPlot() {
        Plot p = new Plot(title:'This bar serves great chips')
        p << new Bars(
            x:(10..70).step(10),
            y:[2,7,1,5,4,6,3], 
            labels:['cat','juice','frog','bog','dog','house','bat']
            )
        p.save('bars.png')
    }

    Random r = new Random()
    
    @Test
    void testHistogramMulti() {
//        List<Double> data = (1..1000).collect { r.nextGaussian() }
        
        Map dataMap = [ (-3) : 3, (-2) : 5, (-1) : 10, 0 : 20, 1 : 10, 2 : 5, 3 : 3, 4 : 1 ]
        Map dataMap2 = [ (-3) : 1, (-2) : 2, (-1) : 5, 0 : 7, 1 : 15, 2 : 12, 3 : 8, 4 : 3,  2: 1 ]
        
        List data = dataMap.collect { [it.key] * it.value }.sum() 
        List data2 = dataMap2.collect { [it.key] * it.value }.sum() 

        println "Data is: " + data
        println "Data2 is: " + data2
        
        Histogram hist = 
            new Histogram(
                data:[data, data2], 
                binCount:8, 
                rangeMax: 3,
                rangeMin: -2,
                names: ['Reference', 'Evaluation'],
                title:'Foos are really frogjible',
                xLabel: 'Frobbles of the Fribfrob',
                yLabel: 'Fringles of the funglefib')
        hist.save('testhist2.png')
    }
 
//    @Test
    void testHistogram() {
//        List<Double> data = (1..1000).collect { r.nextGaussian() }
        
        Map dataMap = [
            (-3) : 3,
            (-2) : 5,
            (-1) : 10,
            0 : 20,
            1 : 10, 
            2 : 5,
            3 : 3,
            4 : 1
        ]
        
        List data = dataMap.collect {
            [it.key] * it.value
        }.sum()
        
        println "Data is: " + data
        
        Histogram hist = 
            new Histogram(
                data:data, 
                binCount:8, 
                title:'Foos are really frogjible',
                xLabel: 'Frobbles of the Fribfrob',
                yLabel: 'Fringles of the funglefib')
        hist.save('testhist.png')
    }
    
    @Test
    void testDensity() {
        def normals = (1..1000).collect { r.nextGaussian() }
        def line = new Density.Line(data:normals, displayName:'TestLine', color: com.twosigma.beakerx.chart.Color.blue)
        assert line.displayName  == 'TestLine'
    }
    
    @Test
    void testDensityBx() {
        def normals = (1..1000).collect { r.nextGaussian() }
        def normals2 = (1..1000).collect { 2 + r.nextGaussian()*1.3 }
        def plot = new Plot()
        def line = new Density.Area(data:normals, 
            displayName:'Test Density', color: com.twosigma.beakerx.chart.Color.blue)
        
        def line2 = new Density.Area(data:normals2,
            displayName:'Higher Test Density', color: com.twosigma.beakerx.chart.Color.green)

        plot  << line
        plot  << line2

        assert line.displayName  == 'Test Density'
        
        plot.save('density.png')
    }
    
    
//    @Test
    void testColoredLines() {
        Plot p = new Plot()
        p << new Line(color: Color.red, x:[1,2,3], y:[4,5,6], width: 1)
        p.save('test.line.png')
    }
}
