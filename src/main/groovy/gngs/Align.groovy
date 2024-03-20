package gngs
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

import org.biojava.nbio.alignment.NeedlemanWunsch;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.alignment.SmithWaterman
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.alignment.template.Aligner;
import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.core.alignment.template.Profile;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.location.template.Location;
 


/**
 * Utility for helping to perform alignments with BioJava
 * 
 * Simple usage of the form:
 * <pre>
 * Align.global('AATTAATTAA','AATTCATAA').profile
 * AATTAATTAA
 * AATTCAT-AA
 * </pre>
 * Or for a more readable base by base comparison, use {@link #profile}:
 * <pre>
 * Align.profile('AATTAATTAA','AATTCATAA')
 * AATTAATTAA
 *     |     
 * AATTCAT-AA
 * </pre>
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
        
    /**
     * Returns a BioJava aligner with profile computed.
     * 
     * @param queryString
     * @param referenceString
     * 
     * @return  aligner with profile computed, use Aligner.profile to access result
     */
    static Aligner global(String queryString, String referenceString) {
        global([:], queryString, referenceString)
    }
    
    /**
     * Returns a BioJava aligner with profile computed.
     * 
     * @param queryString
     * @param referenceString
     * 
     * @return  aligner with profile computed, use Aligner.profile to access result
     */
    static Aligner global(Map params, String queryString, String referenceString) {
        
        // Override defaults with any parameters passed in
        Map actualParams = ALIGNMENT_DEFAULTS + params
        
        DNASequence query = new DNASequence(queryString);
        DNASequence reference = new DNASequence(referenceString);
        
        GapPenalty gapPenalty = new SimpleGapPenalty((short)5,(short)1);
        
        NeedlemanWunsch<DNASequence, NucleotideCompound> aligner = new NeedlemanWunsch<DNASequence, NucleotideCompound>(query, reference, gapPenalty, actualParams.substitutionMatrix);
        
        Profile<DNASequence, NucleotideCompound> profile = aligner.getProfile();
        
        return aligner
    }
    
    /**
     * Returns a formatted string representing the comparison between the two sequences
     * with '-' representing gaps and '|' joining mismatches.
     * 
     * @param queryString
     * @param referenceString
     * @return
     */
    static String profile(String queryString, String referenceString) {
        
        Aligner aligner = Align.global(queryString, referenceString)
        
        
        return aligner.profile
                      .alignedSequences
                      *.asType(List)
                      .transpose()
                      .collect { base1, base2 ->
                         if(base1 == base2) {
                             [base1,' ',base2]
                         }
                         else
                         if(base1.toString() == '-' || base2.toString() == '-') {
                             [base1,' ',base2]
                         }
                         else {
                             [base1, '|', base2]
                         }
                      }.transpose()
                      *.join('')
                      .join('\n')
    }
    
    /**
     * Perform local alignment between the query and reference string
     * 
     * Local alignment means both the query and reference strings are allowed
     * penalty free clipping of their start and end sequence.
     * 
     * @param params
     * @param queryString
     * @param referenceString
     * 
     * @return aligner with profile computed, use Aligner.profile to access result
     */
    static Aligner local(Map params=[:], String queryString, String referenceString) {
        
        // Override defaults with any parameters passed in
        Map actualParams = ALIGNMENT_DEFAULTS + params
        
        DNASequence query = new DNASequence(queryString);
        DNASequence reference = new DNASequence(referenceString);
        
        GapPenalty gapPenalty = new SimpleGapPenalty((short)5,(short)1);
        
        SmithWaterman<DNASequence, NucleotideCompound> aligner = new SmithWaterman<DNASequence, NucleotideCompound>(query, reference, gapPenalty, actualParams.substitutionMatrix);
        
        Profile<DNASequence, NucleotideCompound> profile = aligner.getProfile();
        
        return aligner
    }

}
