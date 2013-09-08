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

import java.util.List;


/**
 * An individual in a family pedigree
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class Subject {
    
    String id
    
    Set<String> phenoTypes
}

/**
 * Models a relationship between two individuals in a pedigree
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class Relationship {
    
    boolean transmits
    
    double transmissionFactor
    
    String type
    
    String from
    
    String to
}

/**
 * Represents a set of relationships between individuals in a family
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class Pedigree {
    
    /**
     * Unique identifier for the family
     */
    String id
    
    /**
     * List of subject ids belonging to the family
     */
    List<String> samples
    
    List<Subject> individuals
    
    List<Integer> phenoTypes
    
    String toString() {
        "$id $samples"
    }
    
    @Lazy
    List<String> affected = { samples.grep {phenoTypes[samples.indexOf(it)] > 0 } }()
   
}
