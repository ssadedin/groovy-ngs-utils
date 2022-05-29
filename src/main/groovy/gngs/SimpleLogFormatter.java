/*
 *  Groovy NGS Utils - Some simple utilites for processing Next Generation Sequencing data.
 *
 *  Copyright (C) 2013 Simon Sadedin, ssadedin<at>gmail.com
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
package gngs;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.logging.*;

/**
 * The default Java log former uses a format that is too verbose, so
 * this class replaces it with something more compact.
 */
public class SimpleLogFormatter extends Formatter {
    
    private static final String lineSep = System.getProperty("line.separator");
    
    private static final int DOT = '.';
    
    /**
     * A Custom format implementation that is designed for brevity.
     */
    public String format(LogRecord record) {
        
        DateFormat format = new SimpleDateFormat("h:mm:ss");
    
        String loggerName = record.getLoggerName();
        if(loggerName == null) {
            loggerName = "root";
        }
        
        int dotIndex = loggerName.lastIndexOf(DOT);
        if(dotIndex>=0)
            loggerName = loggerName.substring(dotIndex+1);
            
        StringBuilder output = new StringBuilder()
            .append(loggerName)
            .append("\t[")
            .append(record.getThreadID())
            .append("]\t")
            .append(record.getLevel()).append("\t|")
            .append(format.format(new Date(record.getMillis())))
            .append(' ')
            .append(record.getMessage()).append(' ')
            .append(lineSep);
            
        if(record.getThrown()!=null) {
            StringWriter w = new StringWriter();
            record.getThrown().printStackTrace(new PrintWriter(w));
            output.append("Exception:\n" + w.toString());
        }    
            
        return output.toString();
    }
}
