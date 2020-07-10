# UTILITY FUNCTIONS ----

# using array notation to allow support for na
not_zero_condition <- function(fieldname) {
  paste0("input['", fieldname, "'] != 0")
}

detail_page_visible <- function() {
  getQueryString()$show == 'detail'
}

create_reactive_data <- function() {
  reactive({
    if (detail_page_visible()) {
      content <- getQueryString()$content
      if (content == 'gene') {
        gene_symbol <- getQueryString()$symbol
      } else if (content == 'pathway') {
        pathway_go <- getQueryString()$go
        pathway_row <- pathways %>%
          filter(go == pathway_go)
        pathway_row$data[[1]]$gene
      } else {
        custom_gene_list <- getQueryString()$custom_gene_list
        c(str_split(custom_gene_list, "\\s*,\\s*", simplify = TRUE))
      }
    }
  })
}

# SHINY MODULES ----

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

#Lucky Gene Link----
# module to display a random interesting gene and navigate to the detail screen for that gene
getLuckyLink <- function(id) {
  ns <- NS(id)
  htmlOutput(ns("get_lucky"), inline = TRUE)
}

surprise <- function(surprise_vec) {
  gene_symbol <- sample(surprise_vec, 1)
  gene_symbol_url <- paste0("?show=detail&content=gene&symbol=", gene_symbol)
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
                        <li>A single gene, such as <a href="?show=detail&content=gene&symbol=TP53">TP53</a> or <a href="?show=detail&content=gene&symbol=BRCA1">BRCA1</a></li>
                        <li>A pathway name, such as <a href="?show=search&query=cholesterol">cholesterol</a>, which will lead you to <a href="?show=detail&content=pathway&go=0006695">Cholesterol Biosynthetic Process</a></li>
                        <li>The Gene Ontology biological process identifier, such as <a href="?show=search&query=1901989">1901989</a>, which will find <a href="?show=detail&content=pathway&go=1901989">Pathway: Positive Regulation Of Cell Cycle Phase Transition (GO:1901989)</a></li>
                        <li>A custom list of genes (separated by commas), such as <a href="?show=search&query=BRCA1,%20BRCA2">BRCA1, BRCA2</a>, which will search <a href="?show=detail&content=custom_gene_list&custom_gene_list=BRCA1,BRCA2">a custom gene list</a></li>
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

# SUMMARY SHINY MODULES ----

detailsSummaryPanel <- function (id) {
  ns <- NS(id)
  uiOutput(ns("detail_summary"))
}

gene_summary_details <- function(gene_summary) {
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
    plotOutput(outputId = "cellanatogram"), 
    hr(), 
    dataTableOutput(outputId = "cellanatogram_table")
  )
}

pathway_summary_details <- function(pathways_row) {
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
    plotOutput(outputId = "cellanatogram"), 
    hr(), 
    dataTableOutput(outputId = "cellanatogram_table")
  )
}

gene_list_summary_details <- function(custom_gene_list) {
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
    plotOutput(outputId = "cellanatogram"),
    hr(),
    dataTableOutput(outputId = "cellanatogram_table")
  )
}

# renders details about gene
gene_summary_ui <- function(gene_symbol) {
  result <- tagList()
  if (gene_symbol != '') {
    gene_summary_row <- gene_summary %>%
      filter(approved_symbol == gene_symbol)
    result <- gene_summary_details(gene_summary_row)
  }
  result
}

pathway_summary_ui <- function(pathway_go) {
  pathway_row <- pathways %>%
    filter(go == pathway_go)
  pathway_summary_details(pathway_row)
}

gene_list_summary_ui <- function(custom_gene_list) {
  # Filter out invalid symbols for when a user edits "custom_gene_list" query parameter
  valid_gene_symbols <- gene_summary %>%
    filter(approved_symbol %in% custom_gene_list) %>%
    pull(approved_symbol)
  gene_list_summary_details(valid_gene_symbols)
}

detailsSummaryPanelServer <- function(id, data) {
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
      
      output$detail_summary <- renderUI({
        if (detail_page_visible()) {
          content <- getQueryString()$content
          if (content == 'gene') {
            gene_summary_ui(data())
          } else if (content == 'pathway') {
            pathway_summary_ui(getQueryString()$go)
          } else {
            custom_gene_list <- getQueryString()$custom_gene_list
            gene_list_summary_ui(str_split(custom_gene_list, "\\s*,\\s*", simplify = TRUE))
          }
        }
      })      
    }
  )
}

# SHINY FUNCTIONS ----
gene_query_result_row <- function(row) {
  gene_summary_row <- row$data
  title <- paste0(gene_summary_row["approved_symbol"], ": ", gene_summary_row["approved_name"])
  list(
    h4(
      tags$strong("Gene:"),
      tags$a(title, href=paste0("?show=detail&content=gene&symbol=", gene_summary_row["approved_symbol"]))
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
      tags$a(title, href=paste0("?show=detail&content=pathway&go=", pathways_row$go))
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
    href <- paste0("?show=detail&content=custom_gene_list&", gene_query_param)
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
