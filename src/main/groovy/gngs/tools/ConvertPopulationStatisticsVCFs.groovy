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
 * brings in significant additional dependencies. To enable parquet support and building of this
 * tool, set ENABLE_PARQUET=true in gradle.properties
 * 
 * @author Simon Sadedin
 */
@Log
class ConvertPopulationStatisticsVCFs extends ToolBase {

    private ParquetWriter writer
    
    Schema schema

    @Override
    public void run() {
        
        if(opts.arguments().isEmpty()) {
            parser.usage()
            System.exit(0)
        }
        
        this.schema = this.createSchema()
        
        writer = AvroParquetWriter.builder(new Path(opts.o)).withSchema(schema).build()
       
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

            record.put("contig", v0.getContig())
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
            o 'Output file', args:1, required: false
        }
    }

}
