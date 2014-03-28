// vim: ts=4:sw=4:expandtab:cindent:
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
 * itself a convenience wrapper around OpenCSV)
 * 
 * @author Simon
 */class TSV {
	
	@Delegate
	CsvIterator parser
	
	TSV(String fileName, List columnNames = null) {
		
		if(!new File(fileName).withReader { r ->
			String firstLine = r.readLine()
			if(firstLine.startsWith('#')) {
				def cols = firstLine.substring(1).trim().split("\t")
				parser = CsvParser.parseCsv(columnNames: cols, separator: '\t', r)
				return true
			}
			return false
		}) {
			new File(fileName).withReader { r ->
				parser = CsvParser.parseCsv(r, separator:'\t')
			}
		}
	}
	
	TSV(Reader reader) {
		String line = reader.readLine()
		parser = CsvParser.parseCsv(reader)
	}	
	
	TSV(Reader reader, List<String> columnNames) {
		parser = CsvParser.parseCsv(reader, columnNames: columnNames, readFirstLine: true, separator: '\t')
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
}
