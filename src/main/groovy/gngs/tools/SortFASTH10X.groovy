/*
 *  Groovy NGS Utils - Some simple utilites for processing Next Generation Sequencing data.
 *
 *  Copyright (C) 2018 Simon Sadedin, ssadedin<at>gmail.com
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
package gngs.tools

import gngs.*
import groovy.transform.CompileStatic
import groovy.util.logging.Log
import java.util.zip.GZIPOutputStream

/**
 * Utility to sort interleaved 10X FASTH format
 *
 * @author Damien Zammit
 */
@Log
class SortFASTH10X extends ToolBase {

    static void main(String [] args) {
        cli('SortFASTH10X < in.fasth', args) {}
    }

    @Override
    public void run() {
        log.info "Sorting FASTH..."
        Writer gzOutput = new BufferedOutputStream(new GZIPOutputStream(System.out), 2048*1024).newWriter()

        try {
            run(gzOutput)
        }
        finally {
            gzOutput.close()
        }
    }

    @CompileStatic
    public run(Writer output) {
        FASTQ.sort10X(output)
        log.info "Completed"
    }
}
