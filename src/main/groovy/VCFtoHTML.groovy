import groovy.xml.StreamingMarkupBuilder

int count = 1

Cli cli = new Cli()
cli.with {
    i 'vcf File', args:1, required:true
    o 'Output HTML file', args:1, required:true
}

opts = cli.parse(args)

def js = [
    "http://ajax.aspnetcdn.com/ajax/jQuery/jquery-2.1.0.min.js",
    "http://cdn.datatables.net/1.10.0/js/jquery.dataTables.min.js",
    "http://ajax.googleapis.com/ajax/libs/jqueryui/1.10.4/jquery-ui.min.js",
    "http://cdnjs.cloudflare.com/ajax/libs/jquery-layout/1.3.0-rc-30.79/jquery.layout.min.js",
    "vcf.js"
]

def css = [
    "http://cdn.datatables.net/1.10.0/css/jquery.dataTables.css",
    "//ajax.googleapis.com/ajax/libs/jqueryui/1.10.4/themes/smoothness/jquery-ui.css",
    "http://cdnjs.cloudflare.com/ajax/libs/jquery-layout/1.3.0-rc-30.79/layout-default.css"
]

new File(opts.o).withWriter { w ->
    w.println """
    <html>
        <head>
    """
    w.println css.collect{"<link rel='stylesheet' href='$it'/>"}.join("\n").stripIndent()
    w.println js.collect{"<script type='text/javascript' src='$it'></script>"}.join("\n").stripIndent()
    
     w.println "<script type='text/javascript'>"
    w.println "var variants = ["
    int i=0;
    Variant last = null;
    VCF.parse(opts.i) { v->
        if(i++>0)
            w.println ","
        w.print(groovy.json.JsonOutput.toJson([v.chr, v.pos, v.ref, v.alt, v.qual, v.info.DP ] + v.header.samples.collect { v.sampleDosage(it) }))
        last = v;
    }
    w.println """];

    var samples = ${groovy.json.JsonOutput.toJson(last.header.samples)};

    var variantTable = null;
    \$(document).ready(function() {
        variantTable = \$.VariantTable('variantTable', samples, variants);
    });
    </script>""";
    
   
    w.println """
    </head>
    <body>
    <p>Loading ...</p>
        <div class="ui-layout-north">
            <h1>VCF File ${new File(opts.i).name}</h1>
            <span id=filterOuter><span id=filters> </span></span>
            <div id=filterHelp style='display:none'> </div>
        </div>
        <div class="ui-layout-center">
            <div id=tableHolder>
            </div>
        </div>
        <div class="ui-layout-south">Report created ${new Date()}</div>
        <div class="ui-layout-east" id="ui-layout-east"></div>
    </body>
    </html>
    """
}



//new File(opts.o).withWriter { w ->
//    w << new StreamingMarkupBuilder(/*, escapeAttributes:false*/).bind {
/*
        html {
          head {
              ${js.collect{"<script type='text/javascript' src='$it'></script>"}.join("\n").stripIndent()}
              ${css.collect{"<link rel='stylesheet' href='$it'/>"}.join("\n").stripIndent()}
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
            
                var variants = [ """; 
                def mkpTmp = mkp;
                VCF.parse(opts.i) { v ->
                      mkpTmp.yieldUnescaped(v.toJson() + ",")
                      return false
                } 
                mkp.yieldUnescaped("];")
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
                  VCF.parse(opts.i) { Variant v ->
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
                      false
                  }
              }
            }
            iframe(id:'igv', name:'igv', src:'about:blank', style:'display:none')
          }
        }
    }
}
*/