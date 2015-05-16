// vim: ts=4:sw=4:expandtab:cindent:

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

/**
 * This class is just an alias to bring the Graxxia TSV class
 * into the default namespace.
 * 
 * @author Simon
 */
class TSV extends graxxia.TSV {

    public TSV(Map options, Reader r) {
        super(options, r);
    }

    public TSV(Map options, String fileName) {
        super(options, fileName);
    }

    public TSV(Reader reader, List<String> columnNames) {
        super(reader, columnNames);
    }

    public TSV(Reader r) {
        super(r);
    }

    public TSV(String fileName) {
        super(fileName);
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
