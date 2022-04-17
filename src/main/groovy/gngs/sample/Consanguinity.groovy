// vim: ts=4:sw=4:expandtab:cindent:
/*
 * Copyright (c) 2016 MCRI, authors
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
package gngs.sample

enum Consanguinity {
	NOT_CONSANGUINEOUS, CONSANGUINEOUS, SUSPECTED, UNKNOWN
    
    private static Map codes = [
            "0" : NOT_CONSANGUINEOUS,
            "1" : CONSANGUINEOUS,
            "2" : SUSPECTED,
            "8" : UNKNOWN,
            "No" : NOT_CONSANGUINEOUS,
            "Yes": CONSANGUINEOUS,
            "Suspected" : SUSPECTED,
            "Unknown" : UNKNOWN
    ]
	
	static Consanguinity decode(String value) {
        value = value?.trim()
        
        // Not strictly to spec: but this is the only non-optional field of many
        // so by itself it forces you to enter many other columns if it is required
        if(!value)
            return UNKNOWN
        
        if(codes.containsKey(value))
            return codes[value]
         
		throw new IllegalArgumentException("Bad consanguinity value [$value] specified")
	}
}

