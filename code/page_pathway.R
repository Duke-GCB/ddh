pathwayPage <- function (id) {
  ns <- NS(id)
  tagList(
    head_tags,
    ddhNavbarPage( 
      tabPanel("Summary",
               div(querySearchInput(ns("search")), style="float: right"),
               pathwaySummaryText(ns("summary"))),
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

pathwayPageServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      data <- reactive({
          pathway_go <- getQueryString()$go
          pathway_row <- pathways %>%
            filter(go == pathway_go)
          pathway_row$data[[1]]$gene
      })
      querySearchServer("search")
      # Home
      pathwaySummaryTextServer("summary", pathway_go = getQueryString()$go)
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