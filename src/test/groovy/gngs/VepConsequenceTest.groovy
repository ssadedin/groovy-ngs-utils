package gngs

import spock.lang.Specification

class VepConsequenceTest extends Specification {

    def 'severityOf() is correct'() {
        expect:
        VepConsequence.severityOf('feature_truncation') > VepConsequence.severityOf('intergenic_variant')
        VepConsequence.values().max { VepConsequence.severityOf(it.term) } == VepConsequence.TRANSCRIPT_ABLATION
        VepConsequence.values().min { VepConsequence.severityOf(it.term) } == VepConsequence.UNKNOWN
    }
}
