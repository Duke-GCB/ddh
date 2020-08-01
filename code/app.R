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

#read data from create_*.R
gene_summary <- readRDS(here::here(app_data_dir, paste0(release, "_gene_summary.Rds")))
pathways <- readRDS(here::here(app_data_dir, paste0(release, "_pathways.Rds")))

#read data from generate_depmap_data.R
achilles <- readRDS(file=here::here(app_data_dir, paste0(release, "_achilles.Rds")))
expression_join <- readRDS(file=here::here(app_data_dir, paste0(release, "_expression_join.Rds")))

#read data from generate_depmap_stats.R
sd_threshold <- readRDS(file = here::here(app_data_dir, paste0(release, "_sd_threshold.Rds")))
achilles_lower <- readRDS(file = here::here(app_data_dir, paste0(release, "_achilles_lower.Rds")))
achilles_upper <- readRDS(file = here::here(app_data_dir, paste0(release, "_achilles_upper.Rds")))
mean_virtual_achilles <- readRDS(file = here::here(app_data_dir, paste0(release, "_mean_virtual_achilles.Rds")))
sd_virtual_achilles <- readRDS(file = here::here(app_data_dir, paste0(release, "_sd_virtual_achilles.Rds")))

#read data from generate_depmap_tables & pathways.R
master_bottom_table <- readRDS(file=here::here(app_data_dir, paste0(release, "_master_bottom_table.Rds")))
master_top_table <- readRDS(file=here::here(app_data_dir, paste0(release, "_master_top_table.Rds")))
master_positive <- readRDS(file=here::here(app_data_dir, paste0(release, "_master_positive.Rds")))
master_negative <- readRDS(file=here::here(app_data_dir, paste0(release, "_master_negative.Rds")))
surprise_genes <- readRDS(file=here::here(app_data_dir, paste0(release, "_surprise_genes.Rds")))
censor_genes <- readRDS(file=here::here(app_data_dir, paste0(release, "_censor_genes.Rds")))

#read data from generate_subcell_data.R
subcell <- readRDS(file=here::here(app_data_dir, paste0(release, "_subcell.Rds")))

#FUNCTIONS-----
#common functions
source(here::here("code", "fun_tables.R"))
source(here::here("code", "fun_plots.R"))
source(here::here("code", "fun_graphs.R"))
source(here::here("code", "fun_reports.R"), local = TRUE)
source(here::here("code", "fun_text.R"))

#SHINY FUNCTIONS-----
source(here::here("code", "shiny_tables.R"), local = TRUE)
source(here::here("code", "shiny_plots.R"), local = TRUE)
source(here::here("code", "shiny_graphs.R"), local = TRUE)
source(here::here("code", "shiny_reports.R"), local = TRUE)
source(here::here("code", "shiny_text.R"), local = TRUE)

### HEAD
head_tags <- tags$head(includeHTML("gtag.html"),includeScript("returnClick.js"))

### universal elements
main_title <- HTML('<a href="." style="color:black;">Data-Driven Hypothesis</a>')
window_title <- "Data-Driven Hypothesis | A Hirschey Lab Resource"

ddhNavbarPage <- function(...) {
  navbarPage(title = main_title, windowTitle = window_title, ...)
}

### list of all pages rendered by this app
page_names <- list(
  home="home",
  search="search",
  gene="gene",
  pathway="pathway",
  gene_list="gene_list"
)

### HOME Modules ----

#SEARCH BOX----
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

#LUCKY GENE----
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

#EXAMPLE SEARCHES----
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
                        <li>A custom list of genes (separated by commas), such as <a href="?show=search&query=BRCA1,%20BRCA2">BRCA1, BRCA2</a>, which will search <a href="?show=gene_list&query_type=custom_gene_list&custom_gene_list=BRCA1,BRCA2">a custom gene list</a></li>
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

#HOME(landing) PAGE----
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
  func <- query_type_to_query_result_row[[row$query_type]]
  func(row)
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

### SEARCH RESULT ROWS ----

gene_query_result_row <- function(row) {
  gene_summary_row <- row$data
  title <- paste0(gene_summary_row["approved_symbol"], ": ", gene_summary_row["approved_name"])
  list(
    h4(
      tags$strong("Gene:"),
      tags$a(title, href=paste0("?show=gene&query_type=gene&symbol=", gene_summary_row["approved_symbol"]))
    ),
    div(tags$strong("Aka:"), gene_summary_row["aka"]),
    div(tags$strong("Entrez ID:"), gene_summary_row["ncbi_gene_id"]),
    hr()
  )
}

pathway_query_result_row <- function(row) {
  pathways_row <- row$data
  gene_symbols <- lapply(pathways_row$data, function(x) { paste(x$gene, collapse=', ') })
  title <- paste0(pathways_row$pathway, " (GO:", pathways_row$go, ")")
  list(
    h4(
      tags$strong("Pathway:"),
      tags$a(title, href=paste0("?show=pathway&query_type=pathway&go=", pathways_row$go))
    ),
    tags$dl(
      tags$dt("Genes"),
      tags$dd(gene_symbols),
    ),
    hr()
  )
}

gene_list_query_result_row <- function(row) {
  gene_summary_rows <- row$data
  title <- row$key
  
  known_gene_symbols <- gene_summary_rows %>% 
    filter(known == TRUE) %>%
    pull(approved_symbol)
  has_known_gene_symbols <- !is_empty(known_gene_symbols)
  
  unknown_gene_symbols <- gene_summary_rows %>% 
    filter(known == FALSE) %>%
    pull(approved_symbol)
  has_unknown_gene_symbols <- !is_empty(unknown_gene_symbols)
  
  known_gene_symbols_tags <- NULL
  if (has_known_gene_symbols) {
    gene_query_param <- paste0("custom_gene_list=", paste(known_gene_symbols, collapse=","))
    href <- paste0("?show=gene_list&query_type=custom_gene_list&", gene_query_param)
    known_gene_symbols_tags <- list(
      tags$h6("Known Gene Symbols"),
      tags$a(paste(known_gene_symbols, collapse=", "), href=href)
    )
  }
  
  unknown_gene_symbols_tags <- NULL
  if (has_unknown_gene_symbols) {
    unknown_gene_symbols_tags <- list(
      tags$h6("Unknown Gene Symbols"),
      tags$div(paste(unknown_gene_symbols, collapse=", "))
    )
  }
  
  list(
    h4(
      tags$strong("Custom Gene List")#,
      #tags$span(title)
    ),
    known_gene_symbols_tags,
    unknown_gene_symbols_tags,
    hr()
  )
}

# specifies how to render the results for a specific query_type
# functions that generate rows in fun_tables.R eg. gene_list_query_results_table()
query_type_to_query_result_row = list(
  gene=gene_query_result_row,
  pathway=pathway_query_result_row,
  gene_list=gene_list_query_result_row
)

# PAGE MODULES-----
source(here::here("code", "page_gene.R"), local = TRUE) ### GENE PAGE ----

# Create output for our router in main UI of Shiny app.
ui <- shinyUI(
  fluidPage(
    uiOutput("pageUI")
  )
)

pages <- list(
  home=homePage(page_names$home),
  search=searchPage(page_names$search),
  gene=genePage(page_names$gene, type = "gene"), #put var here
  pathway=genePage(page_names$pathway, type = "pathway"),
  gene_list=genePage(page_names$gene_list, type = "gene_list")
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
  genePageServer(page_names$gene, type = "gene")
  genePageServer(page_names$pathway, type = "pathway")
  genePageServer(page_names$gene_list, type = "gene_list")
})

shinyApp(ui, server)
