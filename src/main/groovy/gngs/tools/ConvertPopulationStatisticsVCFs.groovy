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

import org.apache.avro.Schema
import org.apache.avro.generic.GenericData
import org.apache.parquet.avro.AvroParquetWriter
import org.apache.parquet.hadoop.ParquetWriter
import org.apache.parquet.hadoop.metadata.CompressionCodecName
import org.apache.hadoop.fs.Path

import groovy.transform.CompileStatic
import groovy.util.logging.Log
import htsjdk.variant.variantcontext.Genotype
import htsjdk.variant.variantcontext.GenotypesContext
import htsjdk.variant.variantcontext.VariantContext

import static gngs.Sex.*

/**
 * Takes a list of population statistics VCFs created by {@link CreatePopulationStatisticsVCF} and
 * convert them to Parquet format, combining them at the same time.
 * <p>
 * NOTE: this class depends on parquet-avro which is NOT bundled with GNGS by default as it 
 * brings in significant additional dependencies.
 * 
 * @author Simon Sadedin
 */
@Log
class ConvertPopulationStatisticsVCFs extends ToolBase {

    private ParquetWriter writer
    
    Schema schema
    
    boolean useNCBIContigs = false
    
    /**
     * Known IDS for GRCh38 contigs
     */
    final static Map<String,String> GRCH38_CONTIGS = [
        "chr1" : "NC_000001.11",
        "chr2" :"NC_000002.12",
        "chr3" : "NC_000003.12",
        "chr4" : "NC_000004.12",
        "chr5" : "NC_000005.10",
        "chr6" : "NC_000006.12",
        "chr7" : "NC_000007.14", 
        "chr8" : "NC_000008.11",
        "chr9" : "NC_000009.12",
        "chr10" : "NC_000010.11",
        "chr11" : "NC_000011.10",
        "chr12" : "NC_000012.12",
        "chr13" : "NC_000013.11",
        "chr14" : "NC_000014.9",
        "chr15" : "NC_000015.10",
        "chr16" : "NC_000016.10",
        "chr17" : "NC_000017.11",
        "chr18" : "NC_000018.10",
        "chr19" : "NC_000019.10",
        "chr20" : "NC_000020.11",
        "chr21" : "NC_000021.9",
        "chr22" : "NC_000022.11",
        "chrX" : "NC_000023.11",
        "chrY" : "NC_000024.10",
    ]

    @Override
    public void run() {
        
        this.useNCBIContigs = opts.ncbi
        
        if(opts.arguments().isEmpty()) {
            parser.usage()
            System.exit(0)
        }
        
        this.schema = this.createSchema()
        
        if(new File(opts.o).exists())
            new File(opts.o).delete()
        
        writer = AvroParquetWriter
            .builder(new Path(opts.o))
            .withSchema(schema)
            .withCompressionCodec(CompressionCodecName.ZSTD)
            .build()
       
        writer.withCloseable {
            VariantContext lastVariant = null
            ProgressCounter progress = new ProgressCounter(withRate: true, withTime:true, extra: {
                "chr: $lastVariant.contig pos: $lastVariant.start"
            })
            VCFSiteWalker locusWalker = new VCFSiteWalker(opts.arguments())
            locusWalker.walk { List<VariantContext> variants ->
                lastVariant = variants[0]
                processSite(variants)
                progress.count()
            }
            progress.end()
        }
    }
    
    Schema createSchema() {
        String    schemaJson = '''
        {
            "type": "record",
            "name": "Variant",
            "fields": [
                {"name": "contig", "type": "string"},
                {"name": "position", "type": "long"},
                {"name": "ref", "type": "string"},
                {"name": "alt", "type": "string"},
                {"name": "hom_count", "type": "int"},
                {"name": "het_count", "type": "int"}
            ]
        }
        '''.stripIndent()
    
        return new Schema.Parser().parse(schemaJson)
    }
    
    /**
     * Process a single site which may have multiple alleles and write out 
     * the parquet record.
     * 
     * @param variants
     */
    @CompileStatic
    void processSite(List<VariantContext> variants) {
        Map<String, List<VariantContext>> alleles = 
            variants.groupBy { VariantContext ctx -> 
                ctx.getAlternateAllele(0).baseString
            }
            
        for(List<VariantContext> alleleVariants : alleles.values()) {

            def record = new GenericData.Record(schema)
            
            VariantContext v0 = alleleVariants.find { it != null }
            
            int hom_count = (int)alleleVariants.sum { VariantContext v ->
                v.getAttributeAsInt('HOM', 0)
            }

            int het_count = (int)alleleVariants.sum { VariantContext v ->
                v.getAttributeAsInt('HET', 0)
            }

            String contig = v0.getContig()
            if(useNCBIContigs) {
                contig = GRCH38_CONTIGS.get(contig, contig)
            }

            record.put("contig", contig)
            record.put("position", v0.start)
            record.put("ref" , v0.getReference().baseString)
            record.put("alt" , v0.getAlternateAllele(0))
            record.put("hom_count" , hom_count)
            record.put("het_count" , het_count)

            writer.write(record)
        }
        
    }
   
    static void main(String [] args) {
        cli('ConvertPopulationStatisticsVCF <vcf1> <vcf2> <vcf3> ...', args) {
            ncbi 'Write NCBI contig identifiers instead of chr values'
            o 'Output file', args:1, required: false
        }
    }
}
