$(document).keyup(function(event) {
    if (event.keyCode == 13) {
      var activeElementId = document.activeElement.id;
      if (activeElementId.endsWith("gene_or_pathway_text")) {
        var searchButtonId =  activeElementId.replace("-gene_or_pathway_text", "-gene_or_pathway_search");
        $("#" + searchButtonId).click();
      }
    }
});