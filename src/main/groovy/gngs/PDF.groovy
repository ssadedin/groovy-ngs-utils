package gngs
import java.awt.Color;
import java.awt.Image;
import java.awt.Toolkit;
import java.awt.image.ImageObserver;
import java.io.FileOutputStream;
import java.util.Date;

import com.lowagie.text.Anchor;
import com.lowagie.text.BadElementException;
//import com.lowagie.text.BaseColor;
import com.lowagie.text.Chapter;
import com.lowagie.text.Document;
import com.lowagie.text.DocumentException;
import com.lowagie.text.Element;
import com.lowagie.text.Font;
//import com.lowagie.text.List;
import com.lowagie.text.ListItem;
import com.lowagie.text.Paragraph;
import com.lowagie.text.Phrase;
import com.lowagie.text.Section;
import com.lowagie.text.pdf.PdfPCell;
import com.lowagie.text.pdf.PdfPTable;
import com.lowagie.text.pdf.PdfWriter;

/**
 * A light weight PDF Builder for creating PDFs to show results.
 * <p>
 * To create a PDF, instantiate a {@link PDF} object and use the {@link #document(String, Closure)}
 * method to initiate the building of a PDF document.
 * <p>
 * <b>Example</b>:
 * <pre>
 * new PDF().document("test.pdf") {
 * 		title("This is My PDF")
 *      br()
 * 		p("This is a PDF document")
 * 		color("RED") { bold {
 * 			p("This paragraph is red and bold")
 *      } }
 *      
 *      table {
 *      	head {
 *      		cells("Name","Age","Height")
 *          }
 *          cells("Simon", 18, 168)
 *          color("BLUE") { cell("Peter") } // Peter will be blue!
 *          cells(28, 92)
 *      }
 *      
 *      // Note: the image will be autoscaled and centered
 *      img("test.jpg")
 * }
 * </pre>
 * <p>
 * <b>Note</b>: This implementation is based on iText, a PDF generation library that
 * changed from MPL to LGPL then to AGPL. This implementation is compatible with the 
 * original MPL licensed version. Some people consider it unwise or unethical to use 
 * that version but this is obviously an individual judgement to be made.
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class PDF {
	
	private static Font defaultFont = new Font(Font.HELVETICA, 12);
		
	List<com.lowagie.text.Element> elementStack = []
	
	Font currentFont = null
	
	List<Font> fontStack = []
	
	PDF() {
		fontStack.push(defaultFont)
	}
	
	void document(File file, Closure c) {
        document(file.absolutePath,c)
    }
    
	void document(String fileName, @DelegatesTo(PDF) Closure c) {
		Document document = new Document();
		PdfWriter.getInstance(document, new FileOutputStream(fileName));
		document.open();
		elementStack.add(new Paragraph())
		
		c.delegate = this
		c()
		
		document.add(elementStack[-1]);
		document.newPage()
		document.close();
	}
	
	void font(Color color, int size, Closure c) {
		Font base = fontStack[-1]
		c.delegate = this
		this.fontStack.push(new Font(base.family, size, base.style, color))
		c()
		this.fontStack.pop()
	}
	
	void fontSize(float size, Closure c) {
		Font base = fontStack[-1]
		c.delegate = this
		this.fontStack.push(new Font(base.family, size, base.style, base.color))
		c()
		this.fontStack.pop()
	}
	
	void color(String color, Closure c) {
		Font base = fontStack[-1]
		c.delegate = this
		Color colorObj = toColor(color)
		this.fontStack.push(new Font(base.family, base.size, base.style, colorObj))
		c()
		this.fontStack.pop()
	}
	
	void p(String text) {
	    def e = elementStack[-1]
        if(e instanceof PdfPCell)
            e.addElement(new Paragraph(text, fontStack[-1]));
        else
    	    elementStack[-1].add(new Paragraph(text, fontStack[-1]));
	}
	
	void bold(Closure c) {
		Font base = fontStack[-1]
		this.fontStack.push(new Font(base.family, base.size, Font.BOLD, base.color))
		c.delegate = this
		c()
		this.fontStack.pop()
	}
	
	void red(Closure c) {
		Font base = fontStack[-1]
		this.fontStack.push(new Font(base.family, base.size, base.style, Color.RED))
		c.delegate = this
		c()
		this.fontStack.pop()
	}
	
	PdfPTable currentTable 
	
	List<Map> tableOpts = []
	
	void table(Map props, Closure c) {
		c.delegate = this
		currentTable = new PdfPTable(props.cols);
		tableOpts.push(props)
		c()
		this.elementStack[-1].add(currentTable)
	}
	
	void table(Closure c) {
		
		c.delegate = this
		c()
		
		// t.setBorderColor(BaseColor.GRAY);
		// t.setPadding(4);
		// t.setSpacing(4);
		// t.setBorderWidth(1);

		this.elementStack[-1].add(currentTable)
	}
	
	List<PdfPCell> headerCells
	
	List<Color> tableBackground = []
	
 	void head(Closure c) {
		 headerCells = []
		 c.delegate = this
		 c()
		 currentTable = new PdfPTable(headerCells.size());
		 currentTable.defaultCell.phrase = new Phrase("",fontStack[-1])
		 headerCells.each { 
			PdfPCell cell = new PdfPCell(new Phrase(String.valueOf(it), fontStack[-1]));
			cell.setHorizontalAlignment(Element.ALIGN_CENTER);
			if(tableBackground) {
				cell.setBackgroundColor(tableBackground[-1])
			}
			else {
				cell.setBackgroundColor(Color.LIGHT_GRAY)
			}
			applyCellProperties(cell)
			currentTable.addCell(cell)
	     }
		 headerCells = null
	 }
     
	 void cell(Object contents) {
		 
		// Handle null values as blanks for convenience
		if(contents == null)
			contents = ""
			
		if(headerCells != null) {
			headerCells << String.valueOf(contents)
		}
		else {
            
			PdfPCell cell = null
            if(contents instanceof Closure) {
    			cell = new PdfPCell()
                this.elementStack.push(cell)
                contents.delegate = this
                contents = contents()
                this.elementStack.pop()
            }
            else {
                cell = new PdfPCell(new Phrase(String.valueOf(contents), fontStack[-1]));
            }
            
			if(tableBackground) {
				cell.setBackgroundColor(tableBackground[-1])
			}
		    applyCellProperties(cell)
			currentTable.addCell(cell);
		}
	 }
	 
	 void applyCellProperties(PdfPCell cell) {
		if(tableOpts) {
			tableOpts[-1].each { key, value ->
				if(cell.hasProperty(key)) {
					cell[key] = value
				}
			}
		}
	}
	 
    /**
     * Set alignment within table cells
     */
	void align(String value, Closure c) {
		Map newOpts = [horizontalAlignment: toAlignment(value)] 
		if(tableOpts) 
			newOpts = tableOpts[-1]  + newOpts
		tableOpts.push(newOpts)
		c.delegate = this
		c()
		tableOpts.pop()
	}	
	
	int toAlignment(String value) {
		if(value == "center")
			Element.ALIGN_CENTER
		else
		if(value == "left")
			Element.ALIGN_LEFT
		else
		if(value == "right")
			Element.ALIGN_RIGHT
	}
	 
	 void cells(Iterable i) {
		 Iterator iter = i.iterator()
		 while(iter.hasNext())
		 	cell(iter.next())
	 }
	 
	 void cells(Object...values) {
		 values.each { cell(it) }
	 }
	 
	 void bg(String color, Closure c) {
		 background(color,c)
	 }
	 
	 /**
	  * Set background color of table cells
	  */
	 void background(String color, Closure c) {
		 tableBackground.push(toColor(color))
		 c.delegate = this
		 c()
		 tableBackground.pop()
	 }
	 
	 void br(int number=1) {
		 for (int i = 0; i < number; i++) {
		   elementStack[-1].add(new Paragraph(" "));
		 }
	 }
	 
	 void title(Object contents) {
		 fontSize(20) { bold{ p(contents) } }
		 br()
	 }
	 
	 void img(String src) {
		 java.awt.Image awtImage = Toolkit.getDefaultToolkit().createImage(src);
		 com.lowagie.text.Image img = com.lowagie.text.Image.getInstance(awtImage, null);
		 
		 int height = -1
		 height = awtImage.getHeight({ java.awt.Image i, int infoflags, int x, int y, int width, int h ->
			 height = h
			 return false
		 } as ImageObserver)
		 while(height<0)
		 	Thread.sleep(50)
			 
//		 println "Height = ${height}"
		 
		 float aspectRatio = (height/ (float)awtImage.getWidth(null))
//		 println "Aspect ratio = $aspectRatio"
		 img.scaleAbsolute((float)300,(float)(aspectRatio * 300))
		 img.absoluteX = 150
		 elementStack[-1].add(img)
	 }
	 
	 Color toColor(String color) {
		color = color.toUpperCase()
		
		// Check for a built in color with the given name
		def field = java.awt.Color.class.fields.find{ it.name == color }
		if(field)
			return field.get(null) 
		
		// Not found - try to decode as hex
	    return Color.decode(color)		
	 }
}
