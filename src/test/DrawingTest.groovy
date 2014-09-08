import static org.junit.Assert.*;

import org.junit.Test;


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
	}

}
