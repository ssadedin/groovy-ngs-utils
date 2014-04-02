import javax.swing.JFrame;

import de.erichseifert.gral.data.DataTable;
import de.erichseifert.gral.plots.XYPlot
import de.erichseifert.gral.plots.lines.DefaultLineRenderer2D;
import de.erichseifert.gral.plots.lines.LineRenderer;
import de.erichseifert.gral.ui.InteractivePanel

class Plot extends JFrame {
	
	Plot() {
		setDefaultCloseOperation(EXIT_ON_CLOSE);
		setSize(800, 600);		
		DataTable data = new DataTable(Double.class, Double.class);
		for (double x = -5.0; x <= 5.0; x+=0.25) {
			double y = 5.0*Math.sin(x);
			data.add(x, y);
		}
		
		XYPlot plot = new XYPlot(data);
		LineRenderer lines = new DefaultLineRenderer2D();
		plot.setLineRenderer(data, lines);
		getContentPane().add(new InteractivePanel(plot));
	}
	
	public static void main(String [] args) {
		Plot plotFrame = new Plot()
		plotFrame.visible = true
		
	}

}
