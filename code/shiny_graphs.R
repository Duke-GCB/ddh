geneNetworkGraph <- function(id) {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(sliderInput(inputId = ns("deg"),
                               label = "Filter \nConnections (<)",
                               value = 2, min = 1, max = 10),
                   sliderInput(inputId = ns("threshold"),
                               label = "n related \nGenes",
                               value =10, min = 10, max = 20),
                   actionButton(inputId = ns("update"), 
                                label = "Update"),
                   width = 3),
      mainPanel(forceNetworkOutput(outputId = ns("graph")),
                width = 9)
    )
  )
}
  
geneNetworkGraphServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      #establish reactive value
      rv <- reactiveValues(degree = 2, 
                           threshold = 10)
      
      #update value upon call
      observeEvent(input$update, {
        rv$degree <- input$deg
        rv$threshold <- input$threshold
      })
      
      output$graph <- renderForceNetwork({
        validate(
          need(data() %in% colnames(achilles), "No data found."))
        withProgress(message = 'Running fancy algorithms', detail = 'Hang tight for 10 seconds', value = 1, {
          make_graph(master_top_table, master_bottom_table, data(), threshold = rv$threshold, deg = rv$degree)
        })
      })
    }
  )
}

