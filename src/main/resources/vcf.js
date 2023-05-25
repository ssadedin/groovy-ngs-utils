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


/**************************** PanelApp API Wrapper ******************************/

class PanelAppApi {
  constructor(api_url = "https://panelapp.agha.umccr.org/api/v1/") {
    this.api_url = new URL(api_url);
    this.panel_infos = [];
  }

  concatted_url(suffix) {
    return new URL(suffix, this.api_url.href).href;
  }

  async get_all_page_items(uri) {
    let data = [];
    // Need to keep iterating the pages until we get a 404
    let currPage = 1;
    while (true) {
      try {
        let response = await $.get(uri);
        data.push(...response.results);
        if(!response.next)
           break
        uri = response.next

      } catch (err) {
        // Exit once we hit the page that can't give us data
        console.log(err)
        break;
      }
      currPage += 1;
    }
    return data;
  }

  async get_panel_infos() {
    if (this.panel_infos.length > 0) {
      return this.panel_infos;
    }
    let uri = this.concatted_url("panels/");
    this.panel_infos = await this.get_all_page_items(uri);
    return this.panel_infos;
  }
  
  get_inheritance(moi) {
    if(moi.includes('MONOALLELIC')) {
        return 'AD'
    }
    else
    if(moi.includes('BIALLELIC')) {
        return 'AR'
    }
    else
    if(moi.includes('X-LINKED')) {
        return 'XL'
    }
  }

  async get_panel_genes(panel_pk) {
    let uri = this.concatted_url(`panels/${panel_pk}/genes/`);
    let genes = await this.get_all_page_items(uri);
    let gene_names = [];
    for (let gene_obj of genes) {
      // Extract just the essential info we need
      gene_names.push(
        {
            gene_symbol : gene_obj.gene_data.gene_symbol,
            confidence_level : gene_obj.confidence_level,
            inheritance: this.get_inheritance(gene_obj.mode_of_inheritance)
        });
    }
    return [gene_names, uri];
  }
}


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


var geneListGenes = {}
var geneLists = {}

async function addGeneList(geneListName) {
    $('#panelAppDialog').dialog('close') 

    console.log("Adding genes for gene list name ", geneListName)
    let gl = api.panel_infos.find(gl => gl.name == geneListName)
    
    geneLists[gl.id] = gl
    let genes_info = await api.get_panel_genes(gl.id)
    for(gene_info of genes_info[0]) {
        var geneListObj = geneListGenes[gene_info.gene_symbol] || {}
        gene_info.gene_list_id = gl.id
        geneListObj[gl.id] = gene_info
        geneListGenes[gene_info.gene_symbol] = geneListObj
    }
    filterTable('variantTable')
    renderFilters()
}

async function selectPanel() {
    window.api = new PanelAppApi()
    
    var gene_lists = await api.get_panel_infos()
    
    console.log("got back gene lists", gene_lists)
    
    $('#panelAppDialog').dialog({
        buttons:[{
                    text: "OK",
                    click: () => {
                        addGeneList()
                    }
                }]
    })
    
    $( "#panelBox" ).autocomplete({
          source: gene_lists.map(gl => gl.name),
          change: (event, gl) =>  { addGeneList($('#panelBox')[0].value) },
          close: (event, gl) =>  { addGeneList($('#panelBox')[0].value) },
    })
}

var sampleSubjects  = {};
var sampleCount = 0;
var affectedSampleIndexes = [];
var familyIndexes = {};
var rowProperties = [];
var highlightLink = null;
var highlightTr = null;
var columns = [];
var visibleColumns = [];
 
var enableGraying = false;

/* Tag Support */
var allTags = [];
var userAnnotations = { 
    tags : {} 
};


// (function ($, window, document, undefined ) {

    var variantTable = null;
    var tableData = []
   
    // Array of strings to 'eval' as filters
    var filters = [];
    var layout = null;
    var familyIndex = null;
    var familyNames = [];
    var myTableId = null;
    
    var indexCounter=0
    
    var TAG_INDEX=indexCounter;
    
    var ID_INDEX = -1
    if(window.showId) {
        ID_INDEX=++indexCounter
    }
    
    var CHR_INDEX=++indexCounter;
    var POS_INDEX=++indexCounter;
    var SIZE_INDEX=++indexCounter;
    var REF_INDEX=++indexCounter;
    var ALT_INDEX=++indexCounter;
    var QUAL_INDEX=++indexCounter;
    var DEPTH_INDEX=++indexCounter;
    var VAF_INDEX=++indexCounter;
    var FAMILIES_INDEX=++indexCounter;
    var GQS_INDEX=++indexCounter;
    var GENE_INDEX=++indexCounter;
    var CONS_INDEX=++indexCounter;
    var MAF_INDEX=++indexCounter;
    var CUSTOMINFOS_INDEX=++indexCounter;
    
    var CUSTOMINFOS_TD_INDEX=CUSTOMINFOS_INDEX
   
    columns = columnNames.map(x =>  { 
        return {
            title: x,
            className:'vcfcol',
            visible: !hiddenColumns.includes(x) 
        }; 
    });
    
    // Redefine gene index to be the first visible column named "gene"
    GENE_INDEX = columns.filter(c => c.visible).findIndex(c => c.title.toLowerCase().endsWith('gene'))
    
    if(customInfos.length>0) {
        CUSTOMINFOS_TD_INDEX=columns.filter(c => c.visible).findIndex(c => c.title == customInfos[0].infoField)
    }
    
    console.log('CUSTOMINFOS_INDEX=' +CUSTOMINFOS_INDEX)
    
    var nonSampleColumnCount = CUSTOMINFOS_INDEX+customInfos.length; // TODO - make this not hard coded!
    var AD_INDEX = -1; // set later because sampleCount is not known until initialised
    
    
    $.VariantTable = function(tableId, samples, variants) {
        
        myTableId = tableId;
        
        sampleSubjects = indexSubjects(pedigrees);
        AD_INDEX = nonSampleColumnCount + sampleCount;

        familyIndexes = indexFamilies(pedigrees);
        for(f in familyIndexes) {
        	familyNames.push(f);
        }
        
        affectedSampleIndexes = findAffectedSamples(samples);
        
    	var tableObj = {};
    
        console.log("init");
    
       
        $('#tableHolder').html('<table cellpadding="0" cellspacing="0" border="0" class="display" id="'+tableId+'"></table>' );
        
        console.log('creating table ...');

        newTable = variants;
        variantTable = $('#'+tableId).DataTable({ 
            "iDisplayLength": 40,
            columns: columns,
            data: variants,
            createdRow: (row,data,dataIndex) => createVariantRow(tableId, row, data, dataIndex)
            });
        
        console.log('done...');
    
        tableData = variantTable.rows().data();
    
        layout = $('body').layout({ 
          applyDefaultStyles: true
        });
        layout.sizePane("north",120);
        
        window.layout = layout;
        
        $('#'+tableId).on('order.dt',  add_display_events);
    
        add_display_events();
    
        // renderFilters();
    
        with($('#filterButtons')) {
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
            
            $button('Preload IGV').click(function() {
                let filteredData = filterTable(tableId)
                let regions = filteredData.map(variant => {
                    let chr = variant[CHR_INDEX]
                    let pos = variant[POS_INDEX]
                    console.log(`Precache ${chr}:${pos}`)
                    return computeIGVRegion(chr, pos, variant[SIZE_INDEX] * 2)
                })
                .join(' ')

                let igvLink = "http://localhost:60151/execute?command=precache+" + encodeURIComponent(regions);
                
                document.getElementById('igvframe').src = igvLink
            })

            $button('PanelApp').click(function() {
                selectPanel()
            })
            
            with($div({id:'panelAppDialog', title: "Select a Panel from PanelApp", style: "display:none"})) {
                $input({id:'panelBox', style:'width:20em'})
            }
            
            
            $span('&nbsp;&nbsp;')

            Object.keys(htmlFilters).forEach( filterName => {
                $button('Filter: ' + filterName).click(() => {
                    filters.push({ expr: htmlFilters[filterName], id: filters.length})
                    renderFilters(tableId)
                    filterTable(tableId)
                })
            })
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
        
        
        initSettings(tableId);
        
        return tableObj;
    };
    
    function saveSettings(prefix) {
        
        if(typeof(prefix) == 'undefined')
                prefix = ''
        
        for(var i=0; i<filters.length; ++i) {
            if(filters[i].id != null)
                filters[i].expr = document.getElementById(filters[i].id).value;
        }
        
        var settingsKey = prefix+":"+location.href;
    
        var variantSettings = {};
        if(typeof(localStorage.variantSettings) != 'undefined')
            variantSettings = JSON.parse(localStorage.variantSettings);
        
        if(!variantSettings[settingsKey])
            variantSettings[settingsKey] = {};
        
        variantSettings[settingsKey].filters = filters;
        variantSettings[settingsKey].userAnnotations = userAnnotations;
        localStorage.variantSettings = JSON.stringify(variantSettings);
    }
    
    function initSettings(tableId, prefix) {
        
        if(typeof(prefix) == 'undefined')
            prefix = ''
        var settingsKey = prefix+":"+location.href;
        
        if(!localStorage.variantSettings) {
             var defaultSettings = { }
             defaultSettings[settingsKey] = { filters: [] }
             console.log("No settings defined: initializing defaults");
             localStorage.variantSettings = JSON.stringify(defaultSettings);
        }
        
       var oldSettings = JSON.parse(localStorage.variantSettings)[settingsKey]; 
       if(typeof(oldSettings) == 'undefined')
           oldSettings = { filters:[]};
       
       loadSettings(tableId, oldSettings);
    }
    
    function loadSettings(tableId, settingsObject) {
        
        console.log("Restoring settings ...");
        
        filters = settingsObject.filters;
        
        if(settingsObject.userAnnotations) {
            userAnnotations = settingsObject.userAnnotations;
            for(var dataIndex in userAnnotations) {
                for(var tag in userAnnotations[dataIndex].tags) {
                    if(allTags.indexOf(tag)<0)
                        allTags.push(tag);
                }
            }
        }
        console.log("User annotations are: " + userAnnotations);
        renderFilters(tableId);
        // updateGeneList();
        filterTable("variantTable");
    }
    
    function initAnnotations(cnvIndex) {
        if(!userAnnotations[cnvIndex])
            userAnnotations[cnvIndex] = { tags: { } };
    }
    
    function addTagToRow(tableId, dataIndex, tr, showFn) {
        var newTag = prompt('Enter the new tag: ','');
        if(!newTag)
            return;
        if(allTags.indexOf(newTag)<0)
            allTags.push(newTag);
                    
        initAnnotations(dataIndex);
        userAnnotations[dataIndex].tags[newTag] = true;
        addRowTags(tr, userAnnotations[dataIndex].tags);
        saveSettings();
        showFn(dataIndex);
        filterTable(tableId);
    }
    
    function addRowTags(row, rowTags) {
        var td0 = row.getElementsByTagName('td')[0];
        
        if(typeof(td0.oldHTML) == 'undefined')
            td0.oldHTML = td0.innerHTML;
        
        var tagDivs = [];
        for(tag in rowTags) {
            tagDivs.push('<span class="rowTagDiv tag'+allTags.indexOf(tag)+'">' + tag + '</span>')
        }
        td0.innerHTML = tagDivs.join(' ') + td0.oldHTML;
    }
    
    function removeRowTags(dataIndex, td,  e) {
        
        e.preventDefault();
        e.stopPropagation();
        
        var tagToRemove = e.target.className.match(/tag[0-9]/)[0];
        
        window.tagToRemove = tagToRemove;
        
        var tagIndex = tagToRemove.replace(/[^0-9]*/,"");
        
        console.log("Remove tag " + tagToRemove + " (index="+tagIndex+") from " + e.target.className + ' from row ' + dataIndex);
        
        delete userAnnotations[dataIndex].tags[allTags[tagIndex]];
        addRowTags(td.parentNode, userAnnotations[dataIndex].tags);
        saveSettings(); 
    }
    
    function indexSubjects(peds) {
        var index = {};
    	for(p in peds) {
            console.log("Indexing " + p);
    		var ped = pedigrees[p];
            for(var i=0; i<ped.length; ++i) {
                var sub = ped[i];
                console.log("Indexing " + sub.id)
            	index[sub.id] = sub;
                sampleCount++;
            }
    	}
        return index;
    }
    
    function findAffectedSamples(samps) {
        var affected = [];
    	for(var i=0; i<samps.length; ++i) {
    		if(sampleSubjects[samples[i]].pheno > 0)
    			affected.push(nonSampleColumnCount+i);
    	}
        return affected;
    }
    
    /**
     * Create a lookup that quickly identifies the columns in the variant table for 
     * proband ("affected"), mother and father of affected in each family of 
     * given hashtable of pedigrees.
     */
    function indexFamilies(peds) {
       var index = {};
       for(ped in peds) {
           var proband = findProband(ped);
           var mother = findByRelationshipToAffected(ped,"MOTHER");
           var father = findByRelationshipToAffected(ped, "FATHER");
    	   index[ped] = {
    		   id : ped,
			   proband : proband ? samples.indexOf(proband.id) + nonSampleColumnCount : null,
        	   mother : mother ? samples.indexOf(mother.id) + nonSampleColumnCount : null,
        	   father : father ? samples.indexOf(father.id) + nonSampleColumnCount : null,
        	   members : $.map(pedigrees[ped], function(sample) { return samples.indexOf(sample.id) + nonSampleColumnCount;})
    	   };
       } 
       return index;
    }
    
    /**
     * Find the individual in family "ped" who has relationship "relationship"
     * to at least one of the affected individuals in the family.
     * 
     * @param 	ped				id of pedigree
     * @param	relationship	type of relationship (MOTHER, FATHER, etc.)
     */
    function findByRelationshipToAffected(ped,relationship) {
        var subjectsFound = $.grep(pedigrees[ped], function(subject,index) {
        	// Look at the 'rel' attribute - is there one with requested relationship type?
        	for(var i=0; i<subject.rel.length; ++i) {
                var rel = subject.rel[i];
        		if(rel.type == relationship) {
        			// is the 'to' affected?
        			if(sampleSubjects[rel.to].pheno>0)
        				return true;
        		}
        	}
            return false;
        });
        if(subjectsFound.length > 0)
    		return subjectsFound[0];
        else
        	return null;
    }
    
    /**
     * Attempt to resolve the 'proband' in the given family. This is
     * taken to be the affected who has the most relationships to other
     * members of the family.
     */
    function findProband(ped) {
        var subjectsFound = $.grep(pedigrees[ped], function(subject,index) {return subject.pheno>0;});
        subjectsFound.sort(function(a,b) { return b.rel.length - a.rel.length;});
        if(subjectsFound.length > 0)
    		return subjectsFound[0];
        else
        	return null;
    }

    function renderFilters(tableId) {
        $('#filters').html('');
        with($('#filters').$span()) {
            console.log("filters = " + filters);
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
        
        $('input').bind('keydown keyup', editing_key_press);
        $('input').each(editing_key_press);
        
        
        $('#activeGeneLists').html(
            Object.values(geneLists).map((gl,i) => { 
                return `<div class="genelist_hd genelist${i}">${gl.name}</div>`
            })
        )
       
        layout.resizeAll();
        layout.sizePane("north",120);
    }
    
    function editing_key_press(e){
        // if(!e.which)editing_restore(this.parentNode);
        var text = $('<span>')
            .html($(this).val())
            .appendTo(this.parentNode);
        var w = text.innerWidth();
        text.remove();
        $(this).width(w+10);
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
    
    function executeFilters(filtersToRun, phase1Data) {
    
        // Update the filters from the fields
        // Filter out rows
        var rowSource;
        var otherRowSource;
        
        var affectedSampleIndex = 0;
        var affectedReferenced = false;
        
        // Family-relative properties
        var familyReferenced = false;
        var familyCountReferenced = false;
        var familyIndex = 0;
        var family = familyIndexes[familyNames[0]];
        var missingFamilyMember = null;
        
        var otherReferenced=false;
        var otherExprs = [];
        var otherIndex = 0;
        var otherData = {
          get gene() { return otherRowSource[GENE_INDEX]; },
          get proband() { familyReferenced=true; if(family.proband==null) missingFamilyMember="PROBAND"; return otherRowSource[family.proband]; },
          get mother() { 
        	  familyReferenced=true; 
        	  otherReferenced=true; 
        	  // console.log("Other Referenced"); 
        	  if(family.mother==null) missingFamilyMember="MOTHER"; 
        	  return otherRowSource[family.mother]; 
          },
          get father() { familyReferenced=true; if(family.father==null) missingFamilyMember="FATHER"; return otherRowSource[family.father]; }
        };
        
        // Other-variant properties
        var data = {
          get chr() { return rowSource[CHR_INDEX]; },
          get pos() { return rowSource[POS_INDEX]; },
          get ref() { return rowSource[REF_INDEX]; },
          get alt() { return rowSource[ALT_INDEX]; },
          get qual() { return rowSource[QUAL_INDEX]; },
          get depth() { return rowSource[DEPTH_INDEX]; },
          get families() { familyCountReferenced=true; return rowSource[FAMILIES_INDEX];},
          get gene() { return rowSource[GENE_INDEX]; },
          get priority() { console.log(genePriorities[rowSource[GENE_INDEX]]); return genePriorities[rowSource[GENE_INDEX]] || 0; },
          get cons() { return rowSource[CONS_INDEX]; },
          get maf() { return rowSource[MAF_INDEX]; },
          get ad() {return rowSource[AD_INDEX]  },
          get vaf() {return rowSource[VAF_INDEX]  },
          
          // Returns the distance from 0.5 (absolute value) for the balance between alleles
          get bal() { return  rowSource[AD_INDEX].map(function(aad) { return ((aad == null) || aad[0] == 0) ? 0 : Math.abs(0.5 - aad[1] /(aad[0] + aad[1])) }) },
          get affected() { affectedReferenced=true; return rowSource[affectedSampleIndexes[affectedSampleIndex]]; },
          get proband() { familyReferenced=true; if(family.proband==null) missingFamilyMember="PROBAND"; return rowSource[family.proband]; },
          get mother() { familyReferenced=true; if(family.mother==null) missingFamilyMember="MOTHER"; return rowSource[family.mother]; },
          get father() { familyReferenced=true; if(family.father==null) missingFamilyMember="FATHER"; return rowSource[family.father]; },
          get other() { otherReferenced=true; return otherData; /*in first phase this is a dummy result, 2nd phase will put in real result*/ }
        };
        
        for(var i=0; i<samples.length;++i) {
            Object.defineProperty(data,samples[i], { get:  partial(function(sampleIndex) { 
                    var dosage = new Number(rowSource[sampleIndex+nonSampleColumnCount]);
                    var ads = rowSource[AD_INDEX][sampleIndex]
                    var gqs = rowSource[GQS_INDEX][sampleIndex]
                    try {
                        dosage.depth = ads.reduce((acc,i) => acc+i);
                        dosage.ad = ads;
                        dosage.vaf = dosage.depth > 0 ? (ads[1] / dosage.depth) : 0
                        dosage.gq = gqs ? gqs[sampleIndex] : 0
                    }
                    catch(e) {
                        console.log(e);
                    }
                    return dosage;
                },i)
            });
        }

        for(var i=0; i<customInfos.length;++i) {
            Object.defineProperty(data,customInfos[i].infoField, { get:  partial(function(customInfoIndex) { 
                    let value = rowSource[customInfoIndex+CUSTOMINFOS_INDEX]
                    return value
                },i)
            });
        }
        
        var newTable = [];
        
        var isPhase2 = typeof(phase1Data)!="undefined";
        var variantData = variants;
        if(isPhase2) {
        	console.log("Executing phase2 on " + phase1Data.length + " variants");
            variantData = phase1Data;
            //return {data:variantData, pendingFilters: []};
        
        }

        rowProperties = new Array(variants.length); // Properties of the rows in the table
        
        var countFilters = 0;
    
        var otherFilterAdded = false;
        for(var index=0; index<variantData.length; ++index) {
            
            rowSource = variantData[index];
            rowSource.index = index+1;
            var thisRowProperties = {};
            
            window.dataval = rowSource;
            window.lastdata = data;
            var pass = true;
            var familyPassed = false; // Set to true if one or more families passes the filter
            var familyCount = 0;
            var familyCountExprs = [];
            with(data) {
                for(var i=0; i<filtersToRun.length; ++i) {
                  var result = false;
                  if(filtersToRun[i].expr == '')
                    result = true;
                  else {
                      otherIndex = (index == 0 ? 1 : 0); // TODO: fails if empty table
                      do { // each family 
                        otherRowSource = variantData[otherIndex];
                        var allOtherIndexesFailed=true;
                        do { // each other variant
                            ++countFilters;
                            var familyResult = true;
                            try { familyResult = eval(filtersToRun[i].expr); } catch(e) { console.log("filter failed: " + e);  familyResult = false; }
                            if(!isPhase2 && otherReferenced) {
                                //console.log("Abort phase1 due to otherReferenced index=" + index);
                                // Can't evaluate this until phase 2 - push it on the list of filters to be
                            	// evaluated in phase2
                                if(!otherFilterAdded) {
                                	otherExprs.push(filtersToRun[i]);
                                    otherFilterAdded = true;
                                }
                                familyPassed = true;
                                familyResult = true;
                                break;
                            }
                            
                            if(!familyResult) {
//                                console.log("Family failed: index=" + index + " otherIndex = " + otherIndex + " family=" + family.id ); // + " proband = " + proband);
                                if(familyCountReferenced) {
                                    // Don't bother evaluating this at all - we'll look at it once all the families are evaluated on other conditions
                                    familyCountExprs.push(filtersToRun[i].expr);
                                	break;
                                }
                                else
                                if(familyReferenced) { // keep searching for a family that passes
                                    // Family has failed - but only for THIS otherIndex. What if another otherIndex lets it pass? Any one will do!
                                }
                                else {
                                    // False result not in context of family (variant level filter) means fail whole filter straight away
                                    // PROBLEM: Here we are failing and it aborts the iteration of otherIndexes too
//                                    console.log("Fail without family contex: abort early");
                                    result = false;
                                	break; 
                                }
                            }
                            else // Pass
                            if(familyReferenced) {
                                console.log("Family passed: index=" + index + " otherIndex = " + otherIndex + " family=" + family.id);
                                allOtherIndexesFailed=false;
                                ++familyCount;
                                familyPassed = true;
                            }
                            
                            ++otherIndex;
                            if(otherIndex == index)
                                ++otherIndex;
                            
                            console.log("Other index = " + otherIndex + " / index = " + index + " / " + variantData.length);
                            otherRowSource = variantData[otherIndex];
                        } while(isPhase2 && otherReferenced && otherIndex<variantData.length);

                        if(allOtherIndexesFailed) {
                            if(enableGraying) {
                                console.log("Gray out family " + family.id + " index = " + index);
                                for(var j=0; j<family.members.length; ++j) {
                                    thisRowProperties[family.members[j]] = { color : 'gray'}; // NOTE: this color is not actually applied, it's hard coded.
                                }
                            }
                        }

                        ++familyIndex;
                        family = familyIndexes[familyNames[familyIndex]];
                        missingFamilyMember = null;
                
                        if(familyReferenced) {
                            for(var k=0;k<familyCountExprs.length;++k) {
                                familyPassed = (familyPassed && eval(familyCountExprs));
                            }
                            familyCountExprs=[];
                            result = familyPassed;
                        }
                        else
                            result = familyResult;
                        
                      otherIndex = (index == 0 ? 1 : 0); // TODO: fails if empty table
                      otherReferenced=false;
                    } while(familyReferenced && familyIndex<familyNames.length);
                    familyIndex=0;
                    family = familyIndexes[familyNames[0]]; 
                    familyReferenced=false;
                  }
                  familyCountReferenced=false;
                  affectedSampleIndex=0;
    
                  if(!result) {
                    pass = false;
                    break;
                  }
                }
            }
            
            if(pass) {
            	newTable.push(variantData[index]);
                if(Object.keys(thisRowProperties).length) {
                	rowProperties[newTable.length-1] = thisRowProperties;
                }
                variantData[index][FAMILIES_INDEX] = familyCount;
            }
            thisRowProperties = {};
        }
        console.log("Filters applied: " + countFilters);
        return { data: newTable, pendingFilters: otherExprs};
    }

    var highlightedRow = null;
    function highlightRow(tr) {
        console.log('addign highlight to ' + tr);
            $('.highlighted').removeClass('highlighted');
        $(tr).addClass('highlighted');
        highlightedRow = tr;
    }

    var newTable = null;
    
    function computeIGVRegion(chr, pos, size) {
        let region = `${chr}:${pos}`;
        let padding = (igvPadding > 0 ? igvPadding : 50)
        if(size>50) {
            if(size<50000) {
                // For intermediate sized events, show the region around the event
                // the value of 150 is chosen here because it roughly represents
                // the length of reads so that the user will be able to fully see 
                // most of the read pairs overlapping the event
                region = chr + ':' + (pos - 150) + '-' + (pos + size + 150);
            }
            else {
                // For very large events, just zoom to 100bp around the pos
                region = chr + ':' + (pos - igvPadding) + '-' + (pos + igvPadding);
            }
        }       
        // For small events, show the user's preferred padding
        else {
            region = chr + ':' + (pos - igvPadding) + '-' + (pos + size + igvPadding);
        }
        return region;
    }
    
    function createVariantRow(tableId, row, data, dataIndex ) {
//        console.log("create row, table id = " + tableId);
        var tds = row.getElementsByTagName('td');
        var igvLink 
        let chr = tds[CHR_INDEX].innerHTML
        let pos = parseInt(tds[POS_INDEX].innerHTML,10)
        let size = 0;
        try {
            size = parseInt(tds[SIZE_INDEX].innerHTML, 10)
        }
        catch(e) {
            // ignore
        }
        
        let region = computeIGVRegion(chr, pos , size)

        if(window.bams && window.bams.length>0) {
            igvLink = "http://localhost:60151/load?locus=" + region;
            igvLink += ('&' + bams.map(bam => {
                return 'file=' + encodeURIComponent(bam)
            }).join('&'));
        }
        else {
            igvLink =  "http://localhost:60151/goto?locus="+region
        }
        
        tds[POS_INDEX].innerHTML = 
            "<a id='variant_"+dataIndex+"_detail' class='igvLink' href='"+igvLink + "'>"+ tds[POS_INDEX].innerHTML + "</a>";
        
        $(tds[POS_INDEX]).find('a').click(function(e) { e.stopPropagation(); $(e.target).addClass('wasVisited'); highlightRow(row); });
        
        let ref = tds[REF_INDEX].innerHTML;
        if(ref.length>25)
            tds[REF_INDEX].innerHTML = ref.substring(0,6) + '(' + (ref.length-12) + 'bp)' + ref.substring(ref.length-6,ref.length-1)
            
        let alt = tds[ALT_INDEX].innerHTML;
        if(alt.length>25)
            tds[ALT_INDEX].innerHTML = alt.substring(0,6) + '(' + (alt.length-12) + 'bp)' + alt.substring(alt.length-6,alt.length-1) 
        
        var tags = userAnnotations[dataIndex] && userAnnotations[dataIndex].tags;
        addRowTags(row, tags);
        
        $(tds[0]).dblclick(partial(removeRowTags,dataIndex, tds[0]));
        $(tds[0]).click(function(evt) {
            console.log('Tag event:', evt)
            if(evt.altKey || evt.metaKey) {
                removeRowTags(dataIndex, tds[0], evt);
            }
            else {
                addTagToRow(tableId, dataIndex, evt.target.parentNode, function(variantIndex) { console.log("tag added for variant " + variantIndex) });
            }
        });

        var ads = data[nonSampleColumnCount+sampleCount].map(function(ad) { return (ad[0] == null) ? "." : ad[0] + "/" + (ad[1]+ad[0]); }).join(", ");
        tds[DEPTH_INDEX].title = ads

        if(rowProperties[dataIndex]) {
            if(enableGraying) {
                for(var i in rowProperties[dataIndex]) {
                    tds[parseInt(i)].style.color = '#eee';
                }
            }
        }
        
/*        
        if(newTable[dataIndex][FAMILIES_INDEX] == 0) {
            var fc = familyNames.reduce(function(prev,curr) {
                return prev + (familyIndexes[curr].members.every(function(sampleIndex) { return newTable[dataIndex][sampleIndex]==0;}) ? 0 : 1);
            },0);
            newTable[dataIndex][FAMILIES_INDEX]=fc;
            tds[FAMILIES_INDEX].innerHTML = fc+"";
        }
 */       
        
        for(var i=0; i<customInfos.length; ++i) {
            if(customInfos[i].type == 'list') { 
                let listValues = data[CUSTOMINFOS_INDEX+i] ? data[CUSTOMINFOS_INDEX+i].split(', ') : []
                if(customInfos[i].rendering == 'label') { 
                    tds[CUSTOMINFOS_TD_INDEX + i].innerHTML = 
                        listValues.map( value => `<div style='display: inline-block;' class='rowTagDiv tag0'>${value}</div>`).join('<br>')
                }
                else {
                    let ol = 
                        '<ol>' + listValues.map(v => `<li>${v}</li>`).join('\n') + '</ol>'
                    tds[CUSTOMINFOS_TD_INDEX + i].innerHTML = ol;
                    tds[CUSTOMINFOS_TD_INDEX + i].style = 'text-align: left;'
                }
            }
            else {
                tds[CUSTOMINFOS_TD_INDEX + i].innerHTML = data[CUSTOMINFOS_INDEX+i]
            }
        }

        var gene = tds[GENE_INDEX].innerHTML;
        let geneCardsLink = `<a href="http://www.genecards.org/cgi-bin/carddisp.pl?gene=${gene}#diseases" target=genecards_win>${gene}</a>`
        let geneListTags = geneListGenes[gene] ? Object.values(geneListGenes[gene]).map((gl,i) => `<span class='genelist${i} genelist'>` + geneLists[gl.gene_list_id].name + '</span>').join('') : ''
        $(tds[GENE_INDEX]).html(geneCardsLink + geneListTags)
        
        let genePriority = genePriorities[gene]
        if(genePriority) {
           $(tds[GENE_INDEX]).addClass(`priority-${genePriority}`) 
        }
    }
    
   
    function add_display_events() {
        // $(".cnvRow").unbind('click').click(function() {
    
        console.log("Add display events");
        
        $("#variantTable tbody tr").unbind('click').click(function() {
            var variantIndex = parseInt($(this).find('a')[0].id.match(/variant_([0-9]*)_detail/)[1],10);
            console.log("Displaying variant " + variantIndex);
    
           $('.highlighted').removeClass('highlighted');
    
            highlightTr = this;
            $(highlightTr).addClass('highlighted');
            
            show_variant_details(variantIndex);
        })
    }
    
    function show_variant_details(variantIndex) {
        console.log("I would show variant " + variantIndex)
    }
    
    function filterTable(tableId) {
        
        updateFilters();
        
        var filterResult = executeFilters(filters);
        
        if(filterResult.pendingFilters.length>0) {
            window.phase1 = filterResult.data;
        	filterResult = executeFilters(filterResult.pendingFilters, filterResult.data);
        }
        
        newTable = filterResult.data;
        
        console.log('Creating table ...');
    
        $('#tableHolder').html('<table id='+tableId+' class=stripe></table>');
        $('#'+tableId).dataTable({     data: newTable, 
                                       columns: columns, 
                                       iDisplayLength: 50,
                                       destroy: true,
                                       createdRow: (row,data,dataIndex) => createVariantRow(tableId, row, data, dataIndex)
                                   });
        add_display_events();
        
        return newTable;
    };
    

    	 
// })( jQuery, window, document );

function open_geneinfo(gene) {
    console.log("Opening gene 2 " + gene)
    $('#southContent').html(`
        <iframe style='width: 100%; height: 100%;' src="http://www.genecards.org/cgi-bin/carddisp.pl?gene=${gene}#diseases"></iframe>
    `)
    layout.sizePane("south",500);
    
    return false
}

