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

package trie

import groovy.transform.CompileStatic

class TrieNode<T> {
    
    static boolean verbose = false
    
    static final int DEFAULT_MAX_VALUES = 50
    
    List<T> values = []
    
    TrieNode() {
    }
    
    @CompileStatic
    TrieNode add(String key, T value) {
       
        if(key.size()==0) {
            values.add(value)
            return this
        }
        
        String childPivot = key[0]
        TrieNode childNode = this.children[childPivot]
        if(!childNode) {
            childNode = new TrieNode()
            this.children[childPivot] = childNode
        }
        childNode.add(key.substring(1),value)
        return this
    }
    
    List<T> getAt(String key) {
        
        if(key.isEmpty())
            return values
        
        String pivot = key[0]
        TrieNode node = children[pivot]
        
        if(!node)
            return []
        
        return node.getAt(key.substring(1))       
    }
    
    public static TrieMismatchCosts DEFAULT_COSTS = new TrieMismatchCosts()
    
    @CompileStatic
    List<TrieQuery> startsWith(String key, TrieQuery query, int maxValues = DEFAULT_MAX_VALUES) {
        
        if(key.isEmpty()) {
            if(verbose)
                println "Query: ${query.path} : $query.cigar on node ${getAllValues(maxValues)} (cost = " + query.cost(DEFAULT_COSTS) + ")"
                       
            query.result = this.getAllValues(maxValues)
            return [query]
        }
        
        String pivot = key[0]
        
        String subKey = key.substring(1)
            
        // Mismatches
        List<TrieQuery<T>> result = []
        for(Map.Entry<String,TrieNode> entry in children) {
            
            String nodeValue = entry.key
            TrieNode node = entry.value
            
            int mismatchValue = 0
            if(nodeValue != pivot) 
                mismatchValue = 1
            
            // Mismatches
            List childResults = []
            int newMismatchValue = query.mismatches - mismatchValue
            if(newMismatchValue>=0) {
                childResults = node.startsWith(subKey, query.mismatch(pivot, newMismatchValue), maxValues)
            }
            
            result.addAll(childResults)
            if(result.size() > maxValues)
                return result.take(maxValues)
                
            // Deletions
            List grandChildResults = null
            if(query.deletions>0) {
                grandChildResults = node.children.iterator().collect { grandChildPivot, TrieNode grandChild -> 
                    int gcMismatchValue = query.mismatches
                    if(grandChildPivot != pivot) 
                        gcMismatchValue -= 1
                    
                    if(gcMismatchValue>=0)
                        return grandChild.startsWith(subKey, query.deletion(pivot).mismatch('',gcMismatchValue), maxValues)
                    else
                       return Collections.EMPTY_LIST
                }
                
                if(result.size() + grandChildResults.size() <= maxValues)
                    result.addAll(grandChildResults)
                else
                    result.addAll(grandChildResults.take(maxValues - result.size()))
            }
        }
        
        if(result.size() > maxValues) {
            return result.take(maxValues)
        }
        
        // Insertions
        List<TrieQuery<T>> insertionResult = Collections.EMPTY_LIST
        if(query.insertions > 0)
             insertionResult = this.startsWith(subKey, query.insertion(pivot), maxValues)
       
        List<TrieQuery<T>> finalResults = result + (List)insertionResult.take(maxValues - result.size())
             
        return (List<TrieQuery<T>>)finalResults.flatten()
    }
    
    @CompileStatic
    List<T> getAllValues(int maxValues) {
        List<T> result = []
        for(TrieNode child in children.values()) {
            List<T> childValues = child.getAllValues(maxValues) 
            for(childValue in childValues) {
                result.add(childValue)
                    if(result.size() >= maxValues)
                        return result
            }
        }
        
        int neededValues = maxValues - result.size()
        if(neededValues > this.values.size())
            result.addAll(this.values)
        else
            result.addAll(this.values.take(neededValues))
        
        result.unique()
        return result
    }
    
    TreeMap<String,TrieNode> children = [:]
}


