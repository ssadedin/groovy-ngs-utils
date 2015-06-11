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

import org.biojava3.alignment.NeedlemanWunsch;
import org.biojava3.alignment.SimpleGapPenalty;
import org.biojava3.alignment.SubstitutionMatrixHelper;
import org.biojava3.alignment.template.Aligner;
import org.biojava3.alignment.template.GapPenalty;
import org.biojava3.alignment.template.Profile;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.location.template.Location;
 
/**
 * Utility for helping to perform alignments with BioJava
 * 
 * @author simon.sadedin@mcri.edu.au
 */
class Align {
    
    // Workaround for Biojava on JDK8
    // see https://stackoverflow.com/a/30656277/4975843
    static Location l = Location.EMPTY
    
    static Map ALIGNMENT_DEFAULTS = [ 
        gapOpenPenalty: 5, 
        gapExtensionPenalty: 1, 
        substitutionMatrix : SubstitutionMatrixHelper.getNuc4_4()
    ]
        
    static Aligner global(String queryString, String referenceString) {
        global([:], queryString, referenceString)
    }
    
    static Aligner global(Map params, String queryString, String referenceString) {
        
        // Override defaults with any parameters passed in
        Map actualParams = ALIGNMENT_DEFAULTS + params
        
        DNASequence query = new DNASequence(queryString);
        DNASequence reference = new DNASequence(referenceString);
        
        GapPenalty gapPenalty = new SimpleGapPenalty((short)5,(short)1);
        
        NeedlemanWunsch<DNASequence, NucleotideCompound> aligner = new NeedlemanWunsch<DNASequence, NucleotideCompound>(query, reference, gapPenalty, actualParams.substitutionMatrix);
//        System.out.println("Score is " +   aligner.getScore());
        
        Profile<DNASequence, NucleotideCompound> profile = aligner.getProfile();
//        System.out.println("Alignment = \n" + profile.toString());
        
        return aligner
    }

}
