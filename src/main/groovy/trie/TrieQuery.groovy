package trie

/*
 *  Groovy NGS Utils - Some simple utilites for processing Next Generation Sequencing data.
 *
 *  Copyright (C) 2017 Simon Sadedin, ssadedin<at>gmail.com
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
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.transform.ToString

class TrieMismatchCosts {
    double match = 0.0d
    
    double mismatch = 1.0d
    
    double insertion = 4.0d
    
    double deletion = 2.0d
}

@ToString(includeNames=true, includePackage=false, excludes="parent")
@CompileStatic
class TrieQuery<T> {
    int mismatches
    int insertions
    int deletions
    int level
    String cigar=""
    String pivot
    
    List<T> result
    
    TrieQuery<T> parent
    
    TrieQuery<T> deletion(String pivot) {
        new TrieQuery(mismatches:mismatches, insertions:insertions, deletions:deletions-1, level:level+1,cigar:cigar+"D", parent:this, pivot: pivot)
    }
    
    TrieQuery<T> insertion(String pivot) {
        new TrieQuery(mismatches:mismatches, insertions:insertions-1, deletions:deletions, level: level+1, cigar: cigar+"I",parent:this, pivot: pivot)
    } 
    
    TrieQuery<T> mismatch(String pivot, int value) {
        new TrieQuery(mismatches:value, insertions:insertions,deletions:deletions,level:level+1, cigar: cigar+(value==mismatches?".":"M"), parent:this, pivot: pivot)
    }
    
    @Memoized
    String getPath() {
        TrieQuery parent = this
        StringBuilder result = new StringBuilder()
        while(parent != null) {
            if(parent.pivot != null)
            result.append(parent.pivot)
            parent = parent.parent
        }
        return result.toString().reverse()
    }
    
    @Memoized
    double cost(TrieMismatchCosts costs) {
        double score = (double)(cigar.iterator().sum { Object cigarValue ->
            switch(cigarValue) {
                case '.':
                    return costs.match
                case 'M':
                    return costs.mismatch
                case 'D':
                    return costs.deletion
                case 'I':
                    return costs.insertion
            }
        } ?: 0.0d)
        return score
    }
}