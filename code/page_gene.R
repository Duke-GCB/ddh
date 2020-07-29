summaryPanel <- function (id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("summary")),
  )
}

gene_summary_details <- function(ns, gene_summary) {
  title <- paste0(gene_summary$approved_symbol, ": ", gene_summary$approved_name)
  tagList(
    h3(title),
    h4("Summary"),
    tags$dl(
      tags$dt("Gene"), tags$dd(gene_summary$approved_symbol),
      tags$dt("Name"), tags$dd(gene_summary$approved_name),
      tags$dt("aka"), tags$dd(gene_summary$aka),
      tags$dt("Entrez ID"), tags$dd(gene_summary$ncbi_gene_id),
      tags$dt("Gene Summary"), tags$dd(gene_summary$entrez_summary)
    ), 
    hr()
    
  )
}

gene_summary_ui <- function(ns, gene_symbol) {
  result <- tagList()
  if (gene_symbol != '') {
    gene_summary_row <- gene_summary %>%
      filter(approved_symbol == gene_symbol)
    result <- gene_summary_details(ns, gene_summary_row)
  }
  result
}

geneSummaryPanelServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
    output$summary <- renderUI({
        gene_summary_ui(session$ns, data())
      })
    }
  )
}

genePage <- function (id) {
  ns <- NS(id)
  tagList(
    head_tags,
    ddhNavbarPage( 
      tabPanel("Summary",
               div(querySearchInput(ns("search")), style="float: right"),
               summaryPanel(ns("summary"))),
      tabPanel(title = "Gene"), #change to navbarMenu when you have a submenu
      tabPanel(title = "Protein"),
      tabPanel(title = "Expression", 
               cellAnatogramPlot(ns("exp")), 
               tags$hr(),
               cellAnatogramTable(ns("exp"))),
      tabPanel(title = "Cells"),
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
                 tabPanel("SPathways",
                          similarPathwaysTable(ns("sim"))),
                 tabPanel("Dissimilar",
                          dissimilarGenesTable(ns("dsim"))),
                 tabPanel("DPathways",
                          dissimilarPathwaysTable(ns("dsim"))),
                 tabPanel("Graph",
                          geneNetworkGraph(ns("graph")))
      ),
      tabPanel("Methods",
               includeHTML(here::here("code","methods.html"))),
      tabPanel("Download Report",
               downloadReportPanel(ns("download"))))
  )
}

genePageServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      data <- reactive({
        if (getQueryString()$show == page_names$gene) {
          getQueryString()$symbol
        }
      })
      querySearchServer("search")
      # Home
      geneSummaryPanelServer("summary", data)
      # Expression
      cellAnatogramPlotServer("exp", data)
      cellAnatogramTableServer("exp", data)
      # Cell Dependencies - Plots
      cellDependenciesPlotServer("dep", data)
      cellBinsPlotServer("dep", data)
      cellDepsLinPlotServer("dep", data)
      # Cell Dependencies - Table
      cellDependenciesTableServer("dep", data)
      # Similar - Genes
      similarGenesTableServer("sim", data)
      # Similar - Pathways
      similarPathwaysTableServer("sim", data)
      # Dissimilar - Genes
      dissimilarGenesTableServer("dsim", data)
      # Dissimilar - Pathways
      dissimilarPathwaysTableServer("dsim", data)
      # Graph
      geneNetworkGraphServer("graph", data)
      # Download
      downloadReportPanelServer("download", data)
    }
  )
}
