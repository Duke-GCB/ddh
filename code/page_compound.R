#this page is for the master compound query, including genes, pathways of genes, and custom gene lists
#it includes a specific function variable 'type' to change shiny module that is called (and assoc. fun) for the description page

compoundSummaryText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("compound_summary_title"))),
    h4("Summary"),
    tags$dl(
      tags$dt("MOA"), tags$dd(textOutput(outputId = ns("compound_summary_moa"))),
      tags$dt("CID"), tags$dd(textOutput(outputId = ns("compound_summary_cid")))
    )
  )
}

compoundSummaryTextServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$compound_summary_title <- renderText({
        data()$id
      })
      output$compound_summary_moa <- renderText({
        summary_compound(prism_names, input = data(), var = "moa")
      })
      output$compound_summary_cid <- renderText({
        summary_compound(prism_names, input = data(), var = "cid")
      })
    })
}

moaSummaryText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("moa_summary_title"))),
    h4("Summary"),
    tags$dl(
      tags$dt("compounds"), tags$dd(textOutput(outputId = ns("moa_summary_compounds"))),
    )
  )
}

moaSummaryTextServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$moa_summary_title <- renderText({
        data()$id
      })
      output$moa_summary_compounds <- renderText({
        paste0(data()$compound, collapse = ", ")
      })
    })
}


compoundListSummaryText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("title"))),
  )
}

compoundListSummaryTextServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$title <- renderText({
        data()$id
      })
    })
}

compoundPage <- function (id, type) {
  ns <- NS(id)

  summary_var <- ""
  if (type == "compound") {
    summary_var <- compoundSummaryText(ns("summary"))
  }
  if (type == "moa") {
    summary_var <- moaSummaryText(ns("summary"))
  }
  if (type == "compound_list") {
    summary_var <- compoundListSummaryText(ns("summary"))
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

compoundPageServer <- function(id, type) {
  moduleServer(
    id,
    function(input, output, session) {
      if(type == "compound") {
        data <- reactive({
          compound_name <- getQueryString()$compound
          list(
            type=type,
            id=compound_name,
            compound=compound_name
          )
        })
      }
      if(type == "moa") {
        data <- reactive({
          moa_str <- getQueryString()$moa
          if (!is.null(moa_str)) {
            prism_rows <- prism_names %>%
              filter(moa == moa_str)
            list(
              type=type,
              id=moa_str,
              compound=prism_rows$name
            )
          }
        })
      }
      if(type == "compound_list") {
        data <- reactive({
          custom_compound_list <- getQueryString()$custom_compound_list
          compounds_list <- c(str_split(custom_compound_list, "\\s*,\\s*", simplify = TRUE))
          list(
            type=type,
            id=custom_compound_list,
            compound=compounds_list
          )
        })
      }
      compoundSummaryTextServer("summary", data)
      moaSummaryTextServer("summary", data)
      compoundListSummaryTextServer("summary", data)
      querySearchServer("search")
      # Download
      downloadReportPanelServer("download", type, data)
    }
  )
}
