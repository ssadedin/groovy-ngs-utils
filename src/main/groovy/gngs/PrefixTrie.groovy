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

import groovy.transform.CompileStatic

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

import trie.*

class PrefixTrie<T> {
    
    Map<String, TrieNode<T>> root = new HashMap<String,TrieNode<T>>()
    
    PrefixTrie() {
    }
    
    PrefixTrie add(String key, T value, int leading=0) {
        
        if(TrieNode.verbose)
            println "Add $key (leading = $leading)"
        
        String pivot = key[0]
        TrieNode node = root[pivot]
        if(!node) {
            root[pivot] = node = new TrieNode()
        }
        if(key.size()>1)
            node.add(key.substring(1), value)
        else
            node.values.add(value)
            
        if(leading>0 && key.size()>1)
            add(key.substring(1), value, leading-1)
            
        return this
    }
    
    List<T> getAt(String key) {
        String pivot = key[0]
        TrieNode node = root[pivot]
        
        if(!node)
            return []
        
        return node.getAt(key.substring(1))
    }
    
    @CompileStatic
    List<TrieQuery> query(String key, int mismatches=0, int insertions=0, int deletions=0, int leading=0, int maxValues=TrieNode.DEFAULT_MAX_VALUES) {
        
        List<TrieQuery> result = []
        if(key.isEmpty()) {
            result.addAll(root.values().collect { TrieNode node ->
                new TrieQuery(mismatches:mismatches, insertions:insertions, deletions:deletions, result: node.getAllValues(maxValues))
            })
        }
        else { 
            
            while(leading>=0 && key.size()>0) {
                String pivot = key[0]
                if(TrieNode.verbose) {
                    println "Query $key"
                }
                
                for(Map.Entry<String,TrieNode> entry in root) {
                    String nodeValue = entry.key
                    TrieNode node = entry.value
                    
                    int mismatchValue = 0
                    if(nodeValue != pivot) 
                        mismatchValue = 1
                    
                    TrieQuery defaultQuery =  new TrieQuery(mismatches:mismatches, insertions:insertions, deletions:deletions)   
                    int newMismatchValue = mismatches - mismatchValue
                    if(newMismatchValue>=0) {
                        List<TrieQuery> mismatchResults = 
                            node.startsWith(key.substring(1), defaultQuery.mismatch(nodeValue, newMismatchValue), maxValues)
                        
                        result.addAll(mismatchResults)
                    }
                        
                    if(result.size() > maxValues) {
                        result = result.take(maxValues)
                        break
                    }
                }
                
                
//                result.addAll((List<TrieQuery>)root.collect { String nodeValue, TrieNode node ->
//                    int mismatchValue = 0
//                    if(nodeValue != pivot) 
//                        mismatchValue = 1
//                    
//                    TrieQuery defaultQuery =  new TrieQuery(mismatches:mismatches, insertions:insertions, deletions:deletions)   
//                    int newMismatchValue = mismatches - mismatchValue
//                    if(newMismatchValue>=0)
//                        return node.startsWith(key.substring(1), defaultQuery.mismatch(nodeValue, newMismatchValue), maxValues)
//                    else
//                        return Collections.EMPTY_LIST
//                }.flatten())
                key = key.substring(1)
                --leading
            }
        }
        
        
        TrieMismatchCosts costs = new TrieMismatchCosts()
        
        return result.flatten()
                     .groupBy { q -> ((TrieQuery)q).result }
                     .collect { node, queries -> ((List<TrieQuery>)queries).min { q -> q.cost(costs) } } // minimum cost query to find the same node
                     .sort { q -> ((TrieQuery)q).cost(costs) }
    }
   
    
    /**
     * Return all values that match the given key, allowing for the specified number
     * of mismatches, insertions, deletions and leading optional bases.
     * 
     * @param key
     * @return  list of values found by searching with these parameters
     */
    @CompileStatic
    List<T> startsWith(String key, int mismatches=0, int insertions=0, int deletions=0, int leading=0) {
        
        List<TrieQuery<T>> result = this.query(key, mismatches, insertions, deletions, leading)
        return result.collect { q -> ((TrieQuery)q).result } 
                     .flatten()
                     .unique()
    } 
}
