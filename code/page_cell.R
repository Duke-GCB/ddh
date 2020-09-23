#this page is for the master cell query, including genes, pathways of genes, and custom gene lists
#it includes a specific function variable 'type' to change shiny module that is called (and assoc. fun) for the description page
cellPage <- function (id, type) {
  ns <- NS(id)

  
  tagList(
    head_tags,
    ddhNavbarPage(
      tabPanel("SUMMARY",
               "test"), #summary variable for alt descriptions
      tabPanel("METHODS",
               includeHTML(here::here("code","methods.html"))),
      # tabPanel("DOWNLOADS",
      #          downloadReportPanel(ns("download"))),
      formContent=querySearchInput(ns("search")))
  )
}

cellPageServer <- function(id, type) {
  moduleServer(
    id,
    function(input, output, session) {

      querySearchServer("search")
      # Download
      downloadReportPanelServer("download", type, data)
    }
  )
}
