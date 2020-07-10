$(document).keyup(function(event) {
    if ($("#search-gene_or_pathway_text").is(":focus") && (event.keyCode == 13)) {
        $("#search-gene_or_pathway_search").click();
    }
});