// vim: ts=4:sw=4:expandtab:cindent:
import groovy.transform.CompileStatic;

import java.io.Reader;

import com.xlson.groovycsv.CsvIterator;
/*
 *  Groovy NGS Utils - Some simple utilites for processing Next Generation Sequencing data.
 *
 *  Copyright (C) 2014 Simon Sadedin, ssadedin<at>gmail.com
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
import com.xlson.groovycsv.CsvParser;

/**
 * A convenience wrapper around Groovy-CSV (which is
 * itself a convenience wrapper around OpenCSV).
 * <p>
 * The main functionality added is auto-conversion of column types.
 * When the first line is read, the data is "sniffed" and the type is 
 * inferred in the order that each column can be parsed as:
 * <li>Double
 * <li>Integer
 * <li>Boolean
 * <li>String (default)
 * <p>
 * The parsed lines are then converted to these data types where possible
 * (where not possible, they just leave the value as Strings).
 * 
 * <p>Example</p>
 * <pre>
 * // File is : foo\tbar\t10\t4.2
 * for(line in new TSV("test.tsv")) {
 *    println line[2] + line[3] // Works as Integer and Float
 * }
 * </pre>
 * Passing column names is simple:
 * <pre>
 * for(line in new TSV("test.tsv", columnNames:['name','something','age','height'])) {
 *   println line.height + line.age
 * }
 * </pre>
 * 
 * @author Simon
 */
class TSV implements Iterable {
	
	CsvIterator parser
    
    Closure reader
    
    Map options
    
    TSV(Map options=[:], String fileName) {
        this.reader =  { new File(fileName).newReader() }
        this.options = options
        if(this.options.containsKey('columnNames') && !this.options.containsKey('readFirstLine'))
            this.options.readFirstLine = true
    }
    
    TSV(Map options=[:], Reader r) {
        this.reader =  { return r }
        this.options = options
    }
     
    CsvIterator newIterator() {
        if(!options.containsKey("separator"))
            options.separator = "\t"
        CsvParser.parseCsv(options, reader())
    }
	
    /*
	TSV(String fileName, List columnNames = null) {
		
        // If the first line starts with # then we treat it as column headers
		if(!new File(fileName).withReader { r ->
			String firstLine = r.readLine()
			if(firstLine.startsWith('#')) {
				def cols = firstLine.substring(1).trim().split("\t")
				parser = CsvParser.parseCsv(columnNames: cols, separator: '\t', readFirstLine: true, r)
				return true
			}
			return false
		}) {
			new File(fileName).withReader { r ->
                if(columnNames)
    				parser = CsvParser.parseCsv(r, columns: columnNames, separator:'\t')
                else
    				parser = CsvParser.parseCsv(r, separator:'\t')
			}
		}
	}
    */
	
//	TSV(Reader reader) {
//		String line = reader.readLine()
//		parser = CsvParser.parseCsv(reader)
//	}	
//	
    
	TSV(Reader reader, List<String> columnNames) {
		parser = CsvParser.parseCsv(reader, columnNames: columnNames, readFirstLine: true, separator: '\t')
	}		
    
    Iterator iterator() {
        
        Iterator i = newIterator()
        
        List columnTypes = options.columnTypes
        
        new Iterator() {
            boolean hasNext() {
                i.hasNext()
            }
            
            Object next() {
                def line = i.next()
                
                if(!columnTypes) {
                    columnTypes = [String] * line.values.size()
                    line.values.eachWithIndex  { v, index ->
                        if(v.isInteger())
                            columnTypes[index] = Integer
                        else
                        if(v.isDouble())
                            columnTypes[index] = Double
                        else
                        if(v in ["true","false"])
                            columnTypes[index] = Boolean
                    }
                }
                
                line.values = TSV.this.convertColumns(line.values, columnTypes)
                
                return line
            }
            
            void remove() {
                throw new UnsupportedOperationException()
            }
        }
    }
    
    @CompileStatic
    List<Object> convertColumns(String [] values, List<Class> columnTypes) {
        List<Object> newValues = values as List
        final int numColumns = columnTypes.size()
        for(int index = 0; index<numColumns; ++index) {
            Class type = columnTypes[index]
            if(values.size()<=index)
                return
            try {
                newValues[index] = values[index].asType(type)
            }
            catch(NumberFormatException e) {
                // Ignore
                // This just leaves the unconverted (string) value in 
                // place
            }
        }
        return newValues
    }
	
	static Iterable parse(Closure c = null) {
		Reader r = new InputStreamReader(System.in)
		if(c != null) {
			for(line in CsvParser.parseCsv(r)) {
				c(line)
			}
		}
		else {
			return CsvParser.parseCsv(r)
		}
	}
    
    /**
     * Read the given file and then write it out to standard output
     * after passing to the given closure, allowing modification
     * @return
     */
    static filter(Closure c) {
        
    }
}

class CSV extends TSV {
    
    CSV(Map options=[:], String fileName) {
        super(options + [separator:','],fileName)
    }
    
    CSV(Map options=[:], Reader r) {
        super(options + [separator:','],r)
    }
 }