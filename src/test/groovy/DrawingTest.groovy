import static org.junit.Assert.*;

import org.junit.Test;

import graxxia.*

class DrawingTest {

	@Test
	public void testDrawing() {
		Drawing d = new Drawing("test.png", 1024, 768, 1000,-10,2000,110)
		d.log = true
		d.color("blue")
		d.line(1005,109, 1950,109)
		d.color("red")
		d.line(1005,100, 1950,100)
		d.color("green")
		d.line(1005,80, 1950,80)
		d.color("orange")
		d.line(1005,-9, 1950,-9)
        
        d.color("black")
        d.line(1000, 10, 1200, 10)
        d.line(1100, -10, 1100, 110)
        
        d.text(1100,10, "hello")
        d.text("world")
        
        d.bar(1400..1700, 50, "Hello World!","Hello Mars!")
	}

}
