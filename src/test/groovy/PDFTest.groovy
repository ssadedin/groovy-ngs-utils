import static org.junit.Assert.*;

import org.junit.Test;

import gngs.PDF


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
					background("light_gray") {
						cell("darney")
					}
					
					background("#eebbbb") {
						cells([1,2,3,4])
					}
                    
                    cell {
                        p("multi")
                        p("line")
                    }
					
				}
			}
			
			color("green") {
				p("This paragraph should be green!")
			}
			
			table(cols:2, padding:20) {
				cells("super","padding")
			}
			
			img("tests/test.jpg")
		}
	}

}
