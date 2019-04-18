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
package gngs

import groovy.transform.CompileStatic;

import java.security.MessageDigest;

/**
 * Simple convenience wrapper that calculates the sha1 hash for a class and returns
 * it formatted in a standard way.
 * 
 * @author Simon Sadedin
 */
class Hash {

    public Hash() {
    }
    
    @CompileStatic
    static String sha1(String value, String encoding=null) {
        def messageDigest = MessageDigest.getInstance("SHA1")
        messageDigest.update(encoding ? value.getBytes(encoding) : value.getBytes() );
        return new BigInteger(1, messageDigest.digest()).toString(16).padLeft( 40, '0' )
    }
}
