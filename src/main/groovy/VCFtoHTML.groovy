import groovy.xml.MarkupBuilder

int count = 1

new MarkupBuilder(escapeAttributes:false).html {
  head {
      script(type:'text/javascript') {
          mkp.yieldUnescaped """
          function highlight(tr) {
              var t = document.getElementById('variantsTable');
              var trs = t.getElementsByTagName('tr');
              for(var i=0; i<trs.length; ++i) {
                  trs[i].style.color = 'black';
              }
              tr.style.color = 'red';
          }
      """
      }
  }
  
  body { 
    table('id': 'variantsTable') {
      thead {
          tr {
              th('#')
              th('Location')
              th('Ref')
              th('Alt')
              th('Qual')
              th('IGV')
          }
      }
      
      tbody {
          VCF.parse(args[0]) { Variant v ->
              tr {
                  td("$count")
                  td("$v.chr:$v.pos")
                  td("$v.ref")
                  td("$v.alt")
                  td("$v.qual")
                  td {
                      a(href:"http://localhost:60151/goto?locus=$v.chr:${v.pos}", onclick:'highlight(this.parentNode.parentNode)', target:'igv', 'igv')
                  }
              }
              ++count
          }
      }
    }
    iframe(id:'igv', name:'igv', src:'about:blank', style:'display:none')
  }
}