#this page is for the master cell query, including genes, pathways of genes, and custom gene lists
#it includes a specific function variable 'type' to change shiny module that is called (and assoc. fun) for the description page


cellSummaryText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("cell_summary_title"))),
  )
}

cellSummaryTextServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$cell_summary_title <- renderText({
        paste0(data()$cell_line, collapse = ", ")
      })
    })
}

lineageSummaryText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("lineage_summary_title"))),
  )
}

lineageSummaryTextServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$lineage_summary_title <- renderText({
        #data()$id
        paste0(data()$cell_line, collapse = ", ")
      })
    })
}

cellPage <- function (id, type) {
  ns <- NS(id)

  summary_var <- paste0(c("TODO", type), collapse=":")
  message(type)
  if (type == "cell") {
    summary_var <- cellSummaryText(ns("summary"))
  }
  if (type == "lineage") {
    summary_var <- lineageSummaryText(ns("summary"))
  }
  
  tagList(
    head_tags,
    ddhNavbarPage(
      tabPanel("SUMMARY",
               summary_var),
      tabPanel("METHODS",
               includeHTML(here::here("code","methods.html"))),
      tabPanel("DOWNLOADS",
                downloadReportPanel(ns("download"))),
      formContent=querySearchInput(ns("search")))
  )
}

cellPageServer <- function(id, type) {
  moduleServer(
    id,
    function(input, output, session) {
      if(type == "cell") {
        data <- reactive({
          cell_line <- getQueryString()$cell_line
          list(
            type=type,
            id=cell_line,
            cell_line=cell_line
          )
        })
      }
      if(type == "lineage") {
        data <- reactive({
          lineage_str <- getQueryString()$lineage
          expression_name_row <- expression_names %>%
            filter(lineage == lineage_str)
          list(
            type=type,
            id=lineage_str,
            cell_line=expression_name_row$cell_line
          )
        })
      }
      cellSummaryTextServer("summary", data)
      lineageSummaryTextServer("summary", data)
      querySearchServer("search")
      # Download
      downloadReportPanelServer("download", type, data)
    }
  )
}
