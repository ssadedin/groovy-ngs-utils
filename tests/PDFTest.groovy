import static org.junit.Assert.*;

import org.junit.Test;


class PDFTest {

	@Test
	public void testCreatePDF() {
		
		
		new PDF().document("tests/test.pdf") {
			
			title("My Fabulous Document")
			
			br()
			
			red {
				p("Hi there, this is simon. This is some really long text that should probably spill over and wrap to make a paragraph. I don't know if PDFBox will do that for me or not! I wonder if it will.")
			}
			p("This text is in the default font.")
			
			fontSize(7) {
				p("tiny text")
				table {
					head {
						cell("hello")
						cell("world")
					}
					
					cells("barney", "charney")
					
					cell("aarney")
					cell("darney")
					
					cells([1,2,3,4])
					
				}
			}
			
			color("green") {
				p("This paragraph should be green!")
			}
			
			img("tests/test.jpg")
		}
	}

}
