import gngs.VCF
import gngs.Variant

class TestUpdate {

    static main(args) {
        
        VCF vcf = VCF.parse(args[0],null) { Variant v ->
            v.update {
                v.info.foo = 'bar'
                v.info.bar = 3
                v.info.fubar = 3.0f
            }
        }
        
        vcf.printHeader()
    }

}
