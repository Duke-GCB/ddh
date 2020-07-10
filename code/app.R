library(shiny)
library(shinyWidgets)
library(tidyverse)
library(plotly)
library(networkD3)
library(corrr)
library(here)
library(lubridate)
library(rmarkdown)
library(markdown)
library(tidygraph)
library(ggraph)
library(viridis)
library(cowplot)
library(plotly)
library(DT)
library(future)
library(promises)
library(gganatogram)

render_report_in_background <- FALSE
if (supportsMulticore()) {
  plan(multicore)
  render_report_in_background <- TRUE
}

#LOAD DATA-----
#read current release information
source(here::here("code", "current_release.R"))

# change data_dir to "tests/data" for faster testing during development
# after running create_test_data.R to populate tests/data from data
data_dir <- "data"

#read data from create_*.R
gene_summary <- readRDS(here::here(data_dir, paste0(release, "_gene_summary.Rds")))
pathways <- readRDS(here::here(data_dir, paste0(release, "_pathways.Rds")))

#read data from generate_depmap_data.R
achilles <- readRDS(file=here::here(data_dir, paste0(release, "_achilles.Rds")))
expression_join <- readRDS(file=here::here(data_dir, paste0(release, "_expression_join.Rds")))

#read data from generate_depmap_stats.R
sd_threshold <- readRDS(file = here::here(data_dir, paste0(release, "_sd_threshold.Rds")))
achilles_lower <- readRDS(file = here::here(data_dir, paste0(release, "_achilles_lower.Rds")))
achilles_upper <- readRDS(file = here::here(data_dir, paste0(release, "_achilles_upper.Rds")))
mean_virtual_achilles <- readRDS(file = here::here(data_dir, paste0(release, "_mean_virtual_achilles.Rds")))
sd_virtual_achilles <- readRDS(file = here::here(data_dir, paste0(release, "_sd_virtual_achilles.Rds")))

#read data from generate_depmap_tables & pathways.R
master_bottom_table <- readRDS(file=here::here(data_dir, paste0(release, "_master_bottom_table.Rds")))
master_top_table <- readRDS(file=here::here(data_dir, paste0(release, "_master_top_table.Rds")))
master_positive <- readRDS(file=here::here(data_dir, paste0(release, "_master_positive.Rds")))
master_negative <- readRDS(file=here::here(data_dir, paste0(release, "_master_negative.Rds")))
surprise_genes <- readRDS(file=here::here(data_dir, paste0(release, "_surprise_genes.Rds")))
censor_genes <- readRDS(file=here::here(data_dir, paste0(release, "_censor_genes.Rds")))

#read data from generate_subcell_data.R
subcell <- readRDS(file=here::here(data_dir, paste0(release, "_subcell.Rds")))

#FUNCTIONS-----
#common functions
source(here::here("code", "fun_tables.R"))
source(here::here("code", "fun_plots.R"))
source(here::here("code", "fun_graphs.R"))
source(here::here("code", "fun_reports.R"))

#SHINY FUNCTIONS-----
source(here::here("code", "shiny_util.R"), local = TRUE)
source(here::here("code", "shiny_plots.R"), local = TRUE)
source(here::here("code", "shiny_tables.R"), local = TRUE)
source(here::here("code", "shiny_graphs.R"), local = TRUE)
source(here::here("code", "shiny_reports.R"), local = TRUE)

### HEAD
head_tags <- tags$head(includeHTML("gtag.html"),includeScript("returnClick.js"))

### universal elements
main_title <- HTML('<a href="." style="color:black;">Data-Driven Hypothesis</a>')
window_title <- "Data-Driven Hypothesis | A Hirschey Lab Resource"

ddhNavbarPage <- function(...) {
  navbarPage(title = main_title, windowTitle = window_title, ...)
}

### HOME (landing) PAGE ----

homePage <- function (id) {
  ns <- NS(id)
  tagList(
    head_tags,
    HTML('<center><img src="hex_ddh.png"></center>'),
    tags$div(
      tags$br(),
      HTML('<center>Data-driven hypothesis is a resource developed by the <a href="http://www.hirscheylab.org" style="color:black;">Hirschey Lab</a> for predicting functional relationships for thousands of genes across the human genome.</center>'), 
      tags$br(),
      tags$br()),
    HTML("<center>"),
    querySearchInput("search"), 
    exampleSearchesLink(ns("examples")), 
    ", ", 
    browsePathwaysLink(ns("pathways")),
    ", or",
    getLuckyLink(ns("lucky")),
    HTML("</center>"),
    exampleSearchesPanel(ns("examples")),
    browsePathwaysPanel(ns("pathways")) 
  )
}

homePageServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      getLuckyServer("lucky")
      exampleSearchesLinkServer("examples")
      exampleSearchesPanelServer("examples")
      browsePathwaysLinkServer("pathways")
      browsePathwaysPanelServer("pathways")
    }
  )
}

### SEARCH PAGE ----
searchPage <- function (id) {
  ns <- NS(id)
  tagList(
    head_tags,
    ddhNavbarPage(),
    div(querySearchInput("search"), style="float: right"),
    h3(textOutput("search_title")),
    div(div(h3("Results", class="panel-title"), class="panel-heading"),
        div(uiOutput(ns("genes_search_result")), class="panel-body"),
        class="bg-info panel panel-default"
    )
  )
}

query_result_row <- function(row) {
  if (row$contents == 'gene') {
    gene_query_result_row(row)
  } else if (row$contents == 'pathway') {
    pathway_query_result_row(row)
  } else {
    gene_list_query_result_row(row)
  }
}

searchPageServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      output$search_title <- renderText({
        query <- getQueryString()
        paste0("Search results for '", query$query, "'")
      })
      output$genes_search_result <- renderUI({
        query <- getQueryString()
        if (grepl(',', query$query)) {
          query_results_table <- gene_list_query_results_table(gene_summary, query$query)
        } else {
          query_results_table <- gene_or_pathway_query_results_table(gene_summary, pathways, query$query)
        }
        if (nrow(query_results_table) > 0) {
          apply(query_results_table, 1, query_result_row)
        }
        else {
          "No results found."
        }
      })      
    }
  )
}

### DETAILS PAGE ----

detailPage <- function (id) {
  ns <- NS(id)
  tagList(
    head_tags,
    ddhNavbarPage( 
       tabPanel("Home",
                div(querySearchInput("search"), style="float: right"),
                detailsSummaryPanel(ns("summary"))),
       navbarMenu(title = "Cell Dependencies",
                  tabPanel("Plots",
                           cellDependenciesPlot(ns("dep")),
                           tags$hr(),
                           cellBinsPlot(ns("dep")),
                           tags$hr(),
                           cellDepsLinPlot(ns("dep"))),
                  tabPanel("Table",
                           cellDependenciesTable(ns("dep")))),
       navbarMenu(title = "Similar",
                  tabPanel("Genes",
                           similarGenesTable(ns("sim"))),
                  tabPanel("Pathways",
                           similarPathwaysTable(ns("sim")))),
       navbarMenu(title = "Dissimilar",
                  tabPanel("Genes",
                           dissimilarGenesTable(ns("dsim"))),
                  tabPanel("Pathways",
                           dissimilarPathwaysTable(ns("dsim")))),
       tabPanel("Graph",
                geneNetworkGraph(ns("graph"))),
       tabPanel("Methods",
                includeHTML(here::here("code","methods.html"))),
       tabPanel("Download Report",
                downloadReportPanel(ns("download"))))
  )
}

detailPageServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      data <- create_reactive_data()
      
      # Home
      detailsSummaryPanelServer("summary", data)
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

# Create output for our router in main UI of Shiny app.
ui <- shinyUI(
  fluidPage(
    uiOutput("pageContent")
  )
)

pages <- list(
  home=homePage("home"),
  search=searchPage("search"),
  detail=detailPage("detail")
)

server <- shinyServer(function(input, output, session) {
  output$pageContent <- renderUI({
    query <- getQueryString()
    show_page <- query$show
    if (is.null(show_page)) {
      show_page <- "home"
    }
    pages[show_page]
  })
  querySearchServer("search")  
  homePageServer("home")
  searchPageServer("search")
  detailPageServer("detail")
})

shinyApp(ui, server)
