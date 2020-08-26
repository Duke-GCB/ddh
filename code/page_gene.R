#this page is for the master gene query, including genes, pathways of genes, and custom gene lists
#it includes a specific function variable 'type' to change shiny module that is called (and assoc. fun) for the description page
genePage <- function (id, type) {
  ns <- NS(id)
  
  #this block is the logic to define the summary_var variable, to display the proper summary module
  if(type == "gene"){
    summary_var <- geneSummaryText(ns("summary"))
    gene_var <- geneText(ns("gene_var"))
  } else if (type == "pathway"){
    summary_var <- pathwaySummaryText(ns("summary"))
    gene_var <- nameText(ns("gene_var"))
  } else if (type == "gene_list") {
    summary_var <- geneListSummaryText(ns("summary")) 
    gene_var <- nameText(ns("gene_var"))
  } else {
    stop("call your summary argument")
  }
  
  tagList(
    head_tags,
    ddhNavbarPage( 
      tabPanel("Summary",
               summary_var), #summary variable for alt descriptions
      tabPanel(title = "Gene", 
               gene_var), #change to navbarMenu when you have a submenu
      tabPanel(title = "Protein"),
      navbarMenu(title = "Expression", 
                 tabPanel("Sub-cell", 
                          cellAnatogramPlot(ns("exp")), 
                          tags$hr(),
                          cellAnatogramTable(ns("exp"))), 
                 tabPanel("Cell line",
                          cellExpressionPlot(ns("cell_exp")),
                          tags$br(),
                          fluidRow(tags$strong(plot_cellexp_title), plot_cellexp_legend,
                                   actionLink(inputId = ns("cell_exp_click"), "View raw expression data below")), #add conditional panel for raw data
                          tags$br(), 
                          conditionalPanel(condition = paste0("input['", ns("cell_exp_click"), "'] != 0"), 
                                           cellExpressionTable(ns("cell_exp"))), 
                          tags$br(),
                 ), 
                 tabPanel("Tissue"), 
                 tabPanel("Coexpression")),
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
               downloadReportPanel(ns("download"))),
      formContent=querySearchInput(ns("search")))
  )
}

genePageServer <- function(id, type) {
  moduleServer(
    id,
    function(input, output, session) {

      #this block is the logic to define the summary_var variable, to display the proper summary module
      if(type == "gene"){
        data <- reactive({
          gene_symbol <- getQueryString()$symbol
          list(
            type=type,
            id=gene_symbol,
            gene_symbols=gene_symbol
          )
        })
        summary_var <- geneSummaryTextServer("summary", data)
        gene_var <- geneTextServer("gene_var", data)
      } else if (type == "pathway"){
        data <- reactive({
          pathway_go <- getQueryString()$go
          pathway_row <- pathways %>%
            filter(go == pathway_go)
          list(
            type=type,
            id=pathway_go,
            gene_symbols=pathway_row$data[[1]]$gene)
        })
        summary_var <- pathwaySummaryTextServer("summary", data)
        gene_var <- nameTextServer("gene_var", data)
      } else if (type == "gene_list") {
        data <- reactive({
          custom_gene_list <- getQueryString()$custom_gene_list
          gene_symbols <- c(str_split(custom_gene_list, "\\s*,\\s*", simplify = TRUE))
          list(
            type=type,
            id=custom_gene_list,
            gene_symbols=gene_symbols
          )
        })
        summary_var <- geneListSummaryTextServer("summary", data)
        gene_var <- nameTextServer("gene_var", data)
      } else {
        stop("fix your summary argument")
      }
      
      querySearchServer("search")
      # Home
      summary_var
      # Gene
      gene_var
      # Protein
      # Expression
      cellAnatogramPlotServer("exp", data)
      cellAnatogramTableServer("exp", data)
      cellExpressionPlotServer("cell_exp", data)
      observeEvent(input$cell_exp_click, { #event to store the 'click'
        })
      cellExpressionTableServer("cell_exp", data)
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
      downloadReportPanelServer("download", type, data)
    }
  )
}
