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
data_dir <- "tests/data"

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

not_zero_condition <- function(fieldname) {
  paste0("input['", fieldname, "'] != 0")
}

### list of all pages rendered by this app
page_names <- list(
  home="home",
  search="search",
  gene="gene",
  pathway="pathway",
  gene_list="gene_list"
)

### HOME (landing) PAGE ----
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
  conditionalPanel(condition = not_zero_condition(ns("example_click")), 
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

#Browse Pathways----
# module that displays a table of pathways when an link is clicked

browsePathwaysLink <- function (id) {
  ns <- NS(id)
  actionLink(inputId = ns("pathway_click"), "browse the pathways")
}

browsePathwaysLinkServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      observeEvent(input$pathway_click, {}) #event to store the 'click'
    }
  )
}

browsePathwaysPanel <- function (id) {
  ns <- NS(id)
  conditionalPanel(condition = not_zero_condition(ns("pathway_click")), 
                   tags$br(),
                   h4("GO Biological Processes"),
                   dataTableOutput(outputId = ns("pathway_table")))
}

browsePathwaysPanelServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      output$pathway_table <- DT::renderDataTable({
        DT::datatable(make_pathway_table(pathways) %>% dplyr::rename(Pathway = pathway, GO = go, Genes = genes), 
                      options = list(pageLength = 10))
      })      
    }
  )
}

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
    querySearchInput(ns("search")), 
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
      querySearchServer("search")
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
    div(querySearchInput(ns("search")), style="float: right"),
    h3(textOutput("search_title")),
    div(div(h3("Results", class="panel-title"), class="panel-heading"),
        div(uiOutput(ns("genes_search_result")), class="panel-body"),
        class="bg-info panel panel-default"
    )
  )
}

query_result_row <- function(row) {
  if (row$query_type == 'gene') {
    gene_query_result_row(row)
  } else if (row$query_type == 'pathway') {
    pathway_query_result_row(row)
  } else {
    gene_list_query_result_row(row)
  }
}

searchPageServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      querySearchServer("search")
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

### GENE PAGE ----

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
    hr(),
    plotOutput(outputId = ns("cellanatogram")), 
    hr(), 
    dataTableOutput(outputId = ns("cellanatogram_table"))
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
      output$cellanatogram <- renderPlot({
        validate(
          need(data() %in% subcell$gene_name, "No subcellular location data for this gene."))
        make_cellanatogram(subcell, data())
      })

      output$cellanatogram_table <- DT::renderDataTable({
        validate(
          need(data() %in% subcell$gene_name, ""))
        DT::datatable(make_cellanatogram_table(subcell, data()), 
                      options = list(pageLength = 10))
      })

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

### PATHWAY PAGE ----

pathway_summary_details <- function(ns, pathways_row) {
  gene_symbols <- lapply(pathways_row$data, function(x) { paste(x$gene, collapse=', ') })
  title <- paste0("Pathway: ", pathways_row$pathway, " (GO:", pathways_row$go, ")")
  list(
    h4(
      tags$strong(title),
    ),
    tags$dl(
      tags$dt("Genes"),
      tags$dd(gene_symbols),
      tags$dt("Pathway Description"), 
      tags$dd(pathways_row$def)
    ),
    hr(),
    plotOutput(outputId = ns("cellanatogram")), 
    hr(), 
    dataTableOutput(outputId = ns("cellanatogram_table"))
  )
}

pathway_summary_ui <- function(ns, pathway_go) {
  pathway_row <- pathways %>%
    filter(go == pathway_go)
  pathway_summary_details(ns, pathway_row)
}

pathwaySummaryPanelServer <- function(id, pathway_go, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$cellanatogram <- renderPlot({
        message(data())
        validate(
          need(data() %in% subcell$gene_name, "No subcellular location data for this gene."))
        make_cellanatogram(subcell, data())
      })
      
      output$cellanatogram_table <- DT::renderDataTable({
        validate(
          need(data() %in% subcell$gene_name, ""))
        DT::datatable(make_cellanatogram_table(subcell, data()), 
                      options = list(pageLength = 10))
      })
      
      output$summary <- renderUI({
        pathway_summary_ui(session$ns, pathway_go)
      })
    }
  )
}

pathwayPage <- function (id) {
  ns <- NS(id)
  tagList(
    head_tags,
    ddhNavbarPage( 
      tabPanel("Summary",
               div(querySearchInput(ns("search")), style="float: right"),
               summaryPanel(ns("summary"))),
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

pathwayPageServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      data <- reactive({
        if (getQueryString()$show == page_names$pathway) {
          pathway_go <- getQueryString()$go
          pathway_row <- pathways %>%
            filter(go == pathway_go)
          pathway_row$data[[1]]$gene
        }
      })
      querySearchServer("search")
      # Home
      pathwaySummaryPanelServer("summary", getQueryString()$go, data)
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

### GENE LIST PAGE ----

gene_list_summary_details <- function(ns, custom_gene_list) {
  gene_symbols <- paste(custom_gene_list, collapse=', ')
  title <- paste0("Custom Gene List") #, gene_symbols
  list(
    h4(
      tags$strong(title),
    ),
    tags$dl(
      tags$dt("Genes"),
      tags$dd(gene_symbols),
    ),
    hr(),
    plotOutput(outputId = ns("cellanatogram")),
    hr(),
    dataTableOutput(outputId = ns("cellanatogram_table"))
  )
}

gene_list_summary_ui <- function(ns, custom_gene_list) {
  # Filter out invalid symbols for when a user edits "custom_gene_list" query parameter
  valid_gene_symbols <- gene_summary %>%
    filter(approved_symbol %in% custom_gene_list) %>%
    pull(approved_symbol)
  gene_list_summary_details(ns, valid_gene_symbols)
}

geneListSummaryPanelServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$cellanatogram <- renderPlot({
        message(data())
        validate(
          need(data() %in% subcell$gene_name, "No subcellular location data for this gene."))
        make_cellanatogram(subcell, data())
      })
      
      output$cellanatogram_table <- DT::renderDataTable({
        validate(
          need(data() %in% subcell$gene_name, ""))
        DT::datatable(make_cellanatogram_table(subcell, data()), 
                      options = list(pageLength = 10))
      })
      
      output$summary <- renderUI({
        custom_gene_list_str <- getQueryString()$custom_gene_list
        custom_gene_list <- str_split(custom_gene_list_str, "\\s*,\\s*", simplify = TRUE)
        gene_list_summary_ui(session$ns, custom_gene_list)
      })
    }
  )
}

geneListPage <- function (id) {
  ns <- NS(id)
  tagList(
    head_tags,
    ddhNavbarPage( 
      tabPanel("Summary",
               div(querySearchInput(ns("search")), style="float: right"),
               summaryPanel(ns("summary"))),
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

geneListPageServer <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      data <- reactive({
        if (getQueryString()$show == page_names$gene_list) {
          custom_gene_list <- getQueryString()$custom_gene_list
          c(str_split(custom_gene_list, "\\s*,\\s*", simplify = TRUE))
        }
      })
      querySearchServer("search")
      # Home
      geneListSummaryPanelServer("summary", data)
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
    uiOutput("pageUI")
  )
)

pages <- list(
  home=homePage(page_names$home),
  search=searchPage(page_names$search),
  gene=genePage(page_names$gene),
  pathway=pathwayPage(page_names$pathway),
  gene_list=geneListPage(page_names$gene_list)
)

server <- shinyServer(function(input, output, session) {
  output$pageUI <- renderUI({
    query <- getQueryString()
    show_page <- query$show
    if (is.null(show_page)) {
      show_page <- page_names$home
    }
    pages[show_page]
  })
  homePageServer(page_names$home)
  searchPageServer(page_names$search)
  genePageServer(page_names$gene)
  pathwayPageServer(page_names$pathway)
  geneListPageServer(page_names$gene_list)
})

shinyApp(ui, server)
