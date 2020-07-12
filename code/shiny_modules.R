# shiny modules that are not graphs, plots, reports or tables.

#Generic modules----
notZeroConditionPanel <- function(fieldname, ...) {
  condition_str <- paste0("input['", fieldname, "'] != 0")
  conditionalPanel(condition = condition_str, ...)  
}

#Search Box----
# module to input search term and navigate to the search screen
querySearchInput <- function(id) {
  ns <- NS(id)
  searchInput(
    inputId = ns("gene_or_pathway"),
    placeholder = "genes, pathways, or GO number", #search
    btnSearch = icon("search")
  )
}

querySearchServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      observeEvent(input$gene_or_pathway_search, { 
        updateQueryString(paste0("?show=search&query=", input$gene_or_pathway), mode="push")
      })
    }
  )
}

#Lucky Gene----
# module to display a random interesting gene and navigate to the detail screen for that gene
getLuckyLink <- function(id) {
  ns <- NS(id)
  htmlOutput(ns("get_lucky"), inline = TRUE)
}

surprise <- function(surprise_vec) {
  gene_symbol <- sample(surprise_vec, 1)
  gene_symbol_url <- paste0("?show=gene&query_type=gene&symbol=", gene_symbol)
  return(gene_symbol_url)
}

getLuckyServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      output$get_lucky <- renderUI({
        tags$a(href = surprise(surprise_genes), "get lucky")
      })
    }
  )
}

#Example Searches----
# module to display a random interesting gene and navigate to the detail screen for that gene

exampleSearchesLink <- function(id) {
  ns <- NS(id)
  actionLink(inputId = ns("example_click"), "See example searches")
}

exampleSearchesLinkServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      observeEvent(input$example_click, {}) # event to store the 'click'
    }
  )
}

exampleSearchesPanel <- function(id) {
  ns <- NS(id)
  notZeroConditionPanel(ns("example_click"), 
                   tags$br(),
                   h4("Examples"),
                   HTML('<h5>Search for</h5>
                        <ul>
                        <li>A single gene, such as <a href="?show=gene&query_type=gene&symbol=TP53">TP53</a> or <a href="?show=gene&query_type=gene&symbol=BRCA1">BRCA1</a></li>
                        <li>A pathway name, such as <a href="?show=search&query=cholesterol">cholesterol</a>, which will lead you to <a href="?show=pathway&query_type=pathway&go=0006695">Cholesterol Biosynthetic Process</a></li>
                        <li>The Gene Ontology biological process identifier, such as <a href="?show=search&query=1901989">1901989</a>, which will find <a href="?show=pathway&query_type=pathway&go=1901989">Pathway: Positive Regulation Of Cell Cycle Phase Transition (GO:1901989)</a></li>
                        <li>A custom list of genes (separated by commas), such as <a href="?show=search&query=BRCA1,%20BRCA2">BRCA1, BRCA2</a>, which will search <a href="?show=pathway&query_type=custom_gene_list&custom_gene_list=BRCA1,BRCA2">a custom gene list</a></li>
                       </ul>')
  )
}

exampleSearchesPanelServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
    }
  )
}


