import javax.swing.JFrame;

import de.erichseifert.gral.data.DataSource;
import de.erichseifert.gral.data.DataTable;
import de.erichseifert.gral.data.EnumeratedData;
import de.erichseifert.gral.data.statistics.Histogram1D;
import de.erichseifert.gral.plots.BarPlot;
import de.erichseifert.gral.plots.XYPlot
import de.erichseifert.gral.plots.lines.DefaultLineRenderer2D;
import de.erichseifert.gral.plots.lines.LineRenderer;
import de.erichseifert.gral.ui.InteractivePanel
import de.erichseifert.gral.util.Orientation;

class Plot extends JFrame {
	
	Plot() {
//		setDefaultCloseOperation(EXIT_ON_CLOSE);
		setSize(800, 600);		
	

	}
    
	public static void main(String [] args) {
		Plot plotFrame = new Plot()
		plotFrame.visible = true
        
//       		DataTable data = new DataTable(Double.class, Double.class);
//		for (double x = -5.0; x <= 5.0; x+=0.25) {
//			double y = 5.0*Math.sin(x);
//			data.add(x, y);
//		}
	 
	}
    
    DataTable toDataTable(Iterable i) {
		DataTable data = new DataTable(Double.class, Double.class);
        i.eachWithIndex { y, index -> data.add((double)index,(double)y) }
        return data
    }
}

class LinePlot extends Plot {
    
    LinePlot(Number...values) {
        this(values as List)
    }
    
    LinePlot(double [] values) {
        this(values as List)
    }    
    
    LinePlot(int [] values) {
        this(values as List)
    }    
     
    LinePlot(Iterable i) {
        DataTable data = toDataTable(i)
		XYPlot plot = new XYPlot(data);
		LineRenderer lines = new DefaultLineRenderer2D();
		plot.setLineRenderer(data, lines);
		getContentPane().add(new InteractivePanel(plot));        
        this.visible = true
    }
    
}

class Histogram extends Plot {
    Histogram(Iterable i) {
        DataTable data = new DataTable(Double.class)
        i.each { data.add((double)it) }
        
        def breaks =  [0,1,2,3,4,5,6,7] as Number[] 
        // int breaks = 10
    
        // Create histogram from data
        Histogram1D histogram = new Histogram1D(data, Orientation.VERTICAL, [breaks] as Number[][]);
    
        // Create a second dimension (x axis) for plotting
        DataSource histogram2d = new EnumeratedData(histogram, 0, (breaks.max() - breaks.min()) / (breaks.size()-1));
        
		LineRenderer lines = new DefaultLineRenderer2D();
        
        BarPlot barPlot = new BarPlot(histogram2d)
		getContentPane().add(new InteractivePanel(barPlot));        
        this.visible = true
    }    
}
