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
import groovyx.gpars.actor.Actor
import groovyx.gpars.actor.Actors

/**
 * Estimates the sex of a sample from one or more VCF file and prints to the console.
 * <p>
 * Usage:
 * <pre>gngs Sex &lt;vcf file&gt;</pre>
 * 
 * See {@link gngs.VCF#guessSex(int)} for details of how sex is determined.
 * 
 * @author Simon Sadedin
 */
@Log
class Sex extends ToolBase {
    
    final static String AMELY_BASES = "CCACTTCCTCCTGCTTGGTCTTGTCTGTTGCTGGCCAAGCTTCCAGATGCAGATCAGGAAGTATGGGGGGCAGGGGCCGCAGGGGGAACATTGGAGGCAGAGGTGGCTGTGGCAGCAGGGGCTGCATGGGTTGCACAGGTGGCTGGGGCTGCATGGGCTGGTGAGGCTGTGGCTGAACAGGCTGGGGCTGGAAGGGCTGCTGGGCAGGCAGAGGGAGGTTTGGCTGATGGTGTTGGGTTGGAGTCATGGATTGCTGGCCAGGAACAGGCATCAGTGCTTGCTGGCGGACCCTGGGCTGCTGAGCTGGCACCACTGGGATGTGGTGATGAGACTGCAGGGTGTGAGTCAGGGGGTGCTGTTGGGACACCACG"
    
    final static MALE_AMELY_READ_THRESHOLD = 100
    
    void run() {
        
        if(!opts.arguments()) {
            parser.usage()
            System.exit(1)
        }
        
        for(String vcf in opts.arguments().grep { it.endsWith('.vcf') }) {
            gngs.Sex sex = new VCF(vcf).guessSex()
            println sex
        }
        
        for(String bamPath in opts.arguments().grep { it.endsWith('.bam') }) {
            if(!opts.t)
                throw new IllegalArgumentException('Please provide a target region to ascertain sex from BAM files')
                
            Regions targetRegions = new BED(opts.t).load()
            SexKaryotyper kt = new SexKaryotyper(new SAM(bamPath), targetRegions)
            kt.run()
            println kt.sex
        }        
        
        for(String fastqPath in opts.arguments().grep { it.endsWith('.fastq.gz') }) {
            gngs.Sex sex = guessFastqSex(fastqPath)
            println sex
        }
    }
    
    @CompileStatic
    gngs.Sex guessFastqSex(String fastqPath) {
        
        // Compute a set of 30 base seeds from the sequence of this chrY exon
        // Note we discard the last entry as it may be less than 30 bases
        final List seeds = AMELY_BASES.toList().collate(30,60)*.join('')[1..-2]
        
        final boolean verbose = opts ? opts['v'] : true
        
        ProgressCounter progress = new ProgressCounter(withRate: true, withTime: true, log: log)
        int count = 0
        Actor counter = Actors.actor {
            loop {
                react { String bases ->
                    if(bases == "stop")
                        terminate()
                    else
                    if(seeds.any { bases.contains(it) })
                        ++count 
                }
            }
        }
        
        FASTQ.eachRead(fastqPath) { FASTQRead r ->
            counter << r.bases
            if(verbose)
                progress.count()
        }
        
        counter << "stop"
        counter.join()
        
        if(verbose) {
            progress.end()
            log.info "Observed $count reads carrying sequence from unique Y chromosome gene"
        }
            
        if(count > 30) {
            return gngs.Sex.MALE
        }
        else
            return gngs.Sex.FEMALE
    }
    
    static main(args) {
        cli("Sex [-t <target region>] <vcf file | bam file | fastq> <vcf file | bam file | fastq> ...", args) {
            t 'Target regions to analyse (required for BAM files)', longOpt: 'target', args:1, required: false
            v 'Print verbose information about how sex is being determined'
        }
    }
}
