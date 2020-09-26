#this page is for the master cell query, including genes, pathways of genes, and custom gene lists
#it includes a specific function variable 'type' to change shiny module that is called (and assoc. fun) for the description page

cellSummaryText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("cell_summary_title"))),
    h4("Summary"),
    tags$dl(
      tags$dt("Lineage"), tags$dd(textOutput(outputId = ns("cell_summary_lineage"))),
      tags$dt("Lineage subtype"), tags$dd(textOutput(outputId = ns("cell_summary_lineage_subtype")))
    )
  )
}

cellSummaryTextServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$cell_summary_title <- renderText({
        data()$id
      })
      output$cell_summary_lineage <- renderText({
        summary_cell(expression_names, input = data(), var = "lineage")
      })
      output$cell_summary_lineage_subtype <- renderText({
        summary_cell(expression_names, input = data(), var = "lineage_subtype")
      })
    })
}

lineageSummaryText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("lineage_summary_title"))),
    h4("Summary"),
    tags$dl(
      tags$dt("Cell lines"), tags$dd(textOutput(outputId = ns("lineage_summary_cell_lines"))),
    )
  )
}

lineageSummaryTextServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$lineage_summary_title <- renderText({
        data()$id
      })
      output$lineage_summary_cell_lines <- renderText({
        paste0(data()$cell_line, collapse = ", ")
      })
    })
}

lineageSubtypeSummaryText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("lineage_subtype_summary_title"))),
    h4("Summary"),
    tags$dl(
      tags$dt("Cells"), tags$dd(textOutput(outputId = ns("lineage_subtype_summary_cell_lines"))),
    )
  )
}

lineageSubtypeSummaryTextServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$lineage_subtype_summary_title <- renderText({
        data()$id
      })
      output$lineage_subtype_summary_cell_lines <- renderText({
        paste0(data()$cell_line, collapse = ", ")
      })
    })
}

cellListSummaryText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("title"))),
  )
}

cellListSummaryTextServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$title <- renderText({
        data()$id
      })
    })
}

cellPage <- function (id, type) {
  ns <- NS(id)

  summary_var <- ""
  if (type == "cell") {
    summary_var <- cellSummaryText(ns("summary"))
  }
  if (type == "lineage") {
    summary_var <- lineageSummaryText(ns("summary"))
  }
  if (type == "lineage_subtype") {
    summary_var <- lineageSubtypeSummaryText(ns("summary"))
  }
  if (type == "cell_list") {
    summary_var <- cellListSummaryText(ns("summary"))
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
          if (!is.null(lineage_str)) {
            expression_name_row <- expression_names %>%
              filter(lineage == lineage_str)
            list(
              type=type,
              id=lineage_str,
              cell_line=expression_name_row$cell_line
            )
          }
        })
      }
      if(type == "lineage_subtype") {
        data <- reactive({
          lineage_subtype_str <- getQueryString()$lineage_subtype
          if (!is.null(lineage_subtype_str)) {
            expression_name_row <- expression_names %>%
              filter(lineage_subtype == lineage_subtype_str)
            list(
              type=type,
              id=lineage_subtype_str,
              cell_line=expression_name_row$cell_line
            )
          }
        })
      }
      if(type == "cell_list") {
        data <- reactive({
          custom_cell_list <- getQueryString()$custom_cell_list
          cell_lines <- c(str_split(custom_cell_list, "\\s*,\\s*", simplify = TRUE))
          list(
            type=type,
            id=custom_cell_list,
            cell_line=cell_lines
          )
        })
      }
      cellSummaryTextServer("summary", data)
      lineageSummaryTextServer("summary", data)
      lineageSubtypeSummaryTextServer("summary", data)
      cellListSummaryTextServer("summary", data)
      querySearchServer("search")
      # Download
      downloadReportPanelServer("download", type, data)
    }
  )
}
