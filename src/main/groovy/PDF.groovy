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

class PDF {
	
	private static Font defaultFont = new Font(Font.HELVETICA, 12);
		
	Paragraph paragraph
	
	Font currentFont = null
	
	List<Font> fontStack = []
	
	PDF() {
		fontStack.push(defaultFont)
	}
	
	void document(String fileName, Closure c) {
		Document document = new Document();
		PdfWriter.getInstance(document, new FileOutputStream(fileName));
		document.open();
		paragraph = new Paragraph();
		
		c.delegate = this
		c()
		
		document.add(paragraph);
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
	
	void p(String text) {
	    paragraph.add(new Paragraph(text, fontStack[-1]));
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
	
	void table(Closure c) {
		
		c.delegate = this
		c()
		

		// t.setBorderColor(BaseColor.GRAY);
		// t.setPadding(4);
		// t.setSpacing(4);
		// t.setBorderWidth(1);

		this.paragraph.add(currentTable)
	}
	
	List<PdfPCell> headerCells
 	void head(Closure c) {
		 headerCells = []
		 c.delegate = this
		 c()
		 currentTable = new PdfPTable(headerCells.size());
		 currentTable.defaultCell.phrase = new Phrase("",fontStack[-1])
		 headerCells.each { 
			PdfPCell cell = new PdfPCell(new Phrase(String.valueOf(it), fontStack[-1]));
			cell.setHorizontalAlignment(Element.ALIGN_CENTER);
			cell.setBackgroundColor(Color.GRAY)
			currentTable.addCell(cell)
	     }
		 headerCells = null
	 }
	 
	 void cell(Object contents) {
		if(headerCells != null) {
			headerCells << String.valueOf(contents)
		}
		else {
			PdfPCell cell = new PdfPCell(new Phrase(String.valueOf(contents), fontStack[-1]));
			currentTable.addCell(cell);
		}
	 }
	 
	 void cells(Iterable i) {
		 Iterator iter = i.iterator()
		 while(iter.hasNext())
		 	cell(iter.next())
	 }
	 
	 void cells(Object...values) {
		 values.each { cell(it) }
	 }
	 
	 void br(int number=1) {
		 for (int i = 0; i < number; i++) {
		   paragraph.add(new Paragraph(" "));
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
			 
		 println "Height = ${height}"
		 
		 float aspectRatio = (height/ (float)awtImage.getWidth(null))
		 println "Aspect ratio = $aspectRatio"
		 img.scaleAbsolute(300,(float)aspectRatio * 300)
		 img.absoluteX = 150
		 paragraph.add(img)
	 }
}
