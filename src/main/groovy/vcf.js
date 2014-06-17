// vim: ts=4:sw=4:cindent:expandtab
/*!
 * jQuery create HTML plugin
 * Original author: @efunneko
 * Licensed under the MIT license
 */

// Simple jQuery plugin to make it a bit cleaner to add
// HTML using jQuery. Instead of doing:
//
// $("body").$("<div/>").append(
//   $("<span/>", {id: 'mySpan'}).append(
//     $("<b/>").html("My bold statement!")));
//
// You can do:
//
// $("body").$div().$span({id: 'mySpan'}).b("My bold statement!");

(function ( $, window, document, undefined ) {

    // Create the defaults once
    var pluginName = 'createHtml',
        defaults = {
        },
        elements = [
            "a", "abbr", "address", "area", "article", "aside", "audio", "b", "base",
            "bdi", "bdo", "blockquote", "body", "br", "button", "canvas", "caption",
            "cite", "code", "col", "colgroup", "command", "data", "datalist", "dd",
            "del", "details", "dfn", "div", "dl", "dt", "em", "embed", "fieldset",
            "figcaption", "figure", "footer", "form", "h1", "h2", "h3", "h4", "h5",
            "h6", "head", "header", "hgroup", "hr", "html", "i", "iframe", "img", "input",
            "ins", "kbd", "keygen", "label", "legend", "li", "link", "main", "map", "mark", "math",
            "menu", "meta", "meter", "nav", "noscript", "object", "ol", "optgroup", "option",
            "output", "p", "param", "pre", "progress", "q", "rp", "rt", "ruby", "s",
            "samp", "script", "section", "select", "small", "source", "span", "strong",
            "style", "sub", "summary", "sup", "svg", "table", "tbody", "td", "textarea",
            "tfoot", "th", "thead", "time", "title", "tr", "track", "u", "ul", "var",
            "video", "wbr"
        ],
        methods = { configure: function(options) { } };

    $[pluginName] = function(method, options) {
        if (methods[method]) {
            methods[method](options);
        }
    };

    // Add all element functions
    $.each(elements, function(i, elName) {
        $.fn["$" + elName] = function(content, attrs) {
            if (typeof(content) == 'object' && typeof(attrs) == 'undefined') {
                attrs = content;
                content = undefined;
            }
            var el = $("<" + elName + ">", attrs);
            if (content && content != "") {
                el.html(content);
            }
            $(this).append(el);
            return el;
        }
    });
    
})( jQuery, window, document );


/**************************** Begin VCF.js ******************************/

function partial(fn) {
    var args = Array.prototype.slice.call(arguments);
    args.shift();
    return function() {
        var new_args = Array.prototype.slice.call(arguments);
        args = args.concat(new_args);
        return fn.apply(window, args);
    };
}

(function ( $, window, document, undefined ) {

    var variantTable = null;
    var tableData = []
    //var columnNames = [];
    var columns = [];
    
    // Array of strings to 'eval' as filters
    var filters = [];
    var layout = null;
    
    $.VariantTable = function(tableId, samples, variants) {
        
    	var tableObj = {};
    
        console.log("init");
    
        columns = columnNames.map(function(x) { return {title: x}; });
        
        $('#tableHolder').html('<table cellpadding="0" cellspacing="0" border="0" class="display" id="'+tableId+'"></table>' );
        
        console.log('creating table ...');
        variantTable = $('#'+tableId).DataTable({ 
            "iDisplayLength": 50,
            columns: columns,
            data: variants
        });
        
        console.log('done...');
    
        tableData = variantTable.rows().data();
    
        layout = $('body').layout({ applyDefaultStyles: true });
    
        $('#'+tableId).on('order.dt',  add_display_events);
    
        add_display_events();
    
        // renderFilters();
    
        with($('#filterOuter')) {
            $button({id:'addFilter'}).$span('Add a filter');
            $button({id:'runFilters'}).$span('Run');
            with($button({id:'clearFilters'})) {
                $span('Clear');
                click(function() { filters=[]; renderFilters(tableId); filterTable(tableId); });
            }
            $button('Share').click(function() {
                alert('Copy the text below to share these filters with somebody else:\n\n'+JSON.stringify(JSON.parse(localStorage.tableFilters)[location.href]));
            });
    
            if(typeof(localStorage.tableFilters) != 'undefined') {
                $button('Import').click(function() {
                    var newFilters = prompt("Paste filters from another session below:");
                    if(newFilters) {
                        filters = JSON.parse(newFilters);
                        renderFilters(tableId);
                        filterTable(tableId);
                    }
                });
            }
        }
    
        $('#addFilter').click(function() {
            filters.push({ expr: '', id: filters.length});
            renderFilters(tableId);
        });
        $('#runFilters').click(function() { filterTable(tableId);});
    
//        var attrs = $.map(variants[0], function(v,k) { return k=='seqnames' ? 'chr' : k });
        // $.each(annotations[9],function(k,v) { attrs.push(k); }); 
        var attrs = columnNames;
    
        $('#addFilter').one("click",function() {
            with($('#filterHelp')) {
                $span("Filter attributes: " + attrs.join(","));
            }
            layout.sizePane("north",150);
            $('#filterHelp').slideDown();
        });
    
        if(localStorage.tableFilters) { 
            var oldFilters = JSON.parse(localStorage.tableFilters)[location.href];
            if(oldFilters) {
                filters = oldFilters;
                renderFilters(tableId);
                filterTable(tableId);
            }
        }
    
        $('h3').each(function() {
            $(this).$button('Filter Out').click(function() {
                var cnvId = this.parentNode.id.replace(/[a-z_]/g,'');
                console.log("filter out cnv " + cnvId);
                filters.push({expr: 'index!="'+cnvId+'"'}); 
                filterTable(tableId);
            });
        });
        return tableObj;
    };

    function renderFilters(tableId) {
        $('#filters').html('');
        with($('#filters').$span()) {
            for(var i=0; i<filters.length; ++i) {
                var f = filters[i];
                if(typeof(f.id) == 'undefined')
                    continue;
                $input({type:'text', id: f.id, value: f.expr}).keyup(function(e) {
                        if(e.keyCode == 13)
                            filterTable(tableId);
                }).focus();
            }
        }
    }
    
    function updateFilters() {
        for(var i=0; i<filters.length; ++i) {
            if(filters[i].id != null)
                filters[i].expr = document.getElementById(filters[i].id).value;
        }
        var allFilters = localStorage.tableFilters ? JSON.parse(localStorage.tableFilters) : {};
        allFilters[location.href] = filters;
        localStorage.tableFilters = JSON.stringify(allFilters);
    }
    
    
    var dataMap = {
        'TRUE' : true,
        'FALSE' : false
    };
    
    function sampleDosage(index, variant) {
    	variant[index + 6];
    }

    function filterTable(tableId) {
    
        updateFilters();
    
        // Update the filters from the fields
        // Filter out rows
        var rowSource;
        
        var data = {
          get chr() { return rowSource[0]; },
          get pos() { return rowSource[1]; },
          get ref() { return rowSource[2]; },
          get alt() { return rowSource[3]; },
          get qual() { return rowSource[4]; },
          get depth() { return rowSource[5]; },
          get gene() { return rowSource[6]; },
          get cons() { return rowSource[7]; },
          get maf() { return rowSource[8]; },
        };
        
        for(var i=0; i<samples.length;++i) {
            Object.defineProperty(data,samples[i], { get:  partial(function(sampleIndex) { 
            		return rowSource[sampleIndex+9];
                },i)
            });
        }
        
        var newTable = [];
    
        for(var index=0; index<variants.length; ++index) {
            
            rowSource = variants[index];
//            if(typeof(rowSource.chr) == 'undefined')
//                rowSource.chr = rowSource.seqnames;
            
            rowSource.index = index+1;
            
            window.dataval = rowSource;
            window.lastdata = data;
   
            var indexp = index+1;
//            $('#cnv_'+indexp+'_detail').show();
//            $('#cnv_'+indexp+'_img').show();
    
            var pass = true;
            with(data) {
                for(var i=0; i<filters.length; ++i) {
                  var result;
                  if(filters[i].expr == '')
                    result = true;
                  else {
                    try { result = eval(filters[i].expr); }
                    catch(e) { result = false; }
                  }
    
                  if(!result) {
//                    $('#cnv_'+indexp+'_detail').hide();
//                    $('#cnv_'+indexp+'_img').hide();
                    pass = false;
//                    console.log("fail");
                    break;
                  }
                }
            }
            if(pass)
            	newTable.push(variants[index]);
        }
        
        console.log('Creating table ...');
    
        $('#tableHolder').html('<table id='+tableId+' class=stripe></table>');
        $('#'+tableId).dataTable({ data: newTable, 
                                   columns: columns, 
                                   iDisplayLength: 50,
                                   destroy: true,
                                   });
        add_display_events();
    };
    
    function add_display_events() {
        
    }
    	 
})( jQuery, window, document );
