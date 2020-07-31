#this page is for the master gene query, including genes, pathways of genes, and custom gene lists
#it includes a specific function variable 'summary_text_var' to change shiny module that is called (and assoc. fun) for the description page
genePage <- function (id, summary_text_var) {
  ns <- NS(id)
  
  #this block is the logic to define the summary_var variable, to display the proper summary module
  if(summary_text_var == "gene_summary"){
    summary_var <- geneSummaryText(ns("summary"))
  } else if (summary_text_var == "pathway_summary"){
    summary_var <- pathwaySummaryText(ns("summary"))
  } else if (summary_text_var == "gene_list_summary") {
    summary_var <- geneListSummaryText(ns("summary")) 
  } else {
    stop("call your summary argument")
  }
  
  tagList(
    head_tags,
    ddhNavbarPage( 
      tabPanel("Summary",
               div(querySearchInput(ns("search")), style="float: right"),
               summary_var), #summary variable for alt descriptions
      tabPanel(title = "Gene", 
               geneText(ns("gene_desc"))), #change to navbarMenu when you have a submenu
      tabPanel(title = "Protein"),
      navbarMenu(title = "Expression", 
                 tabPanel("Sub-cell", 
                          cellAnatogramPlot(ns("exp")), 
                          tags$hr(),
                          cellAnatogramTable(ns("exp"))), 
                 tabPanel("Cell line"), 
                 tabPanel("Tissue")),
      navbarMenu(title = "Dependencies",
                 tabPanel("Plots",
                          cellDependenciesPlot(ns("dep")),
                          tags$hr(),
                          cellBinsPlot(ns("dep")),
                          tags$hr(),
                          cellDepsLinPlot(ns("dep"))),
                 tabPanel("Table",
                          cellDependenciesTable(ns("dep"))),
                 tabPanel("Similar",
                          similarGenesTable(ns("sim"))),
                 tabPanel("Dissimilar",
                          dissimilarGenesTable(ns("dsim"))),
                 tabPanel("Graph",
                          geneNetworkGraph(ns("graph")))
      ),
      tabPanel("Methods",
               includeHTML(here::here("code","methods.html"))),
      tabPanel("Download Report",
               downloadReportPanel(ns("download"))))
  )
}

genePageServer <- function(id, summary_text_var) {
  moduleServer(
    id,
    function(input, output, session) {
      #this block is the logic to define the summary_var variable, to display the proper summary module
      if(summary_text_var == "gene_summary"){
        data <- reactive({if (getQueryString()$show == page_names$gene) {getQueryString()$symbol}})
        summary_var <- geneSummaryTextServer("summary", data)
      } else if (summary_text_var == "pathway_summary"){
        data <- reactive({
          pathway_go <- getQueryString()$go
          pathway_row <- pathways %>%
            filter(go == pathway_go)
          pathway_row$data[[1]]$gene
        })
        summary_var <- pathwaySummaryTextServer("summary", pathway_go = getQueryString()$go)
      } else if (summary_text_var == "gene_list_summary") {
        data <- reactive({
          custom_gene_list <- getQueryString()$custom_gene_list
          c(str_split(custom_gene_list, "\\s*,\\s*", simplify = TRUE))
        })
        summary_var <- geneListSummaryTextServer("summary", data)
      } else {
        stop("fix your summary argument")
      }
      
      querySearchServer("search")
      # Home
      summary_var
      # Gene
      geneTextServer("gene_desc", data)
      # Protein
      # Expression
      cellAnatogramPlotServer("exp", data)
      cellAnatogramTableServer("exp", data)
      # Dependencies
      cellDependenciesPlotServer("dep", data)
      cellBinsPlotServer("dep", data)
      cellDepsLinPlotServer("dep", data)
      cellDependenciesTableServer("dep", data)
      # Similar - Genes
      similarGenesTableServer("sim", data)
      # Dissimilar - Genes
      dissimilarGenesTableServer("dsim", data)
      # Graph
      geneNetworkGraphServer("graph", data)
      # Download
      downloadReportPanelServer("download", data)
    }
  )
}
