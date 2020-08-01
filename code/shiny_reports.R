downloadReportPanel <- function(id) {
  ns <- NS(id)
  tagList(
    h2("Report Generator"),    
                     "To generate a report, click on the button below",
                     br(),
                     downloadButton(outputId = ns("report"), label = "Download report")
  )
}

downloadReportPanelServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {   
      output$report <- downloadHandler(
        # create pdf report
        filename = function() {
          type = getQueryString()$content
          if(type == "gene"){
            paste0(data(), "_ddh.pdf")
          } else if (type == "pathway") {
            go <- getQueryString()$go
            paste0("go_", go, "_ddh.pdf")
          } else {
            paste0("custom_", data(), "_ddh.pdf")
          }
        },
        content = function(file) {
          gene_symbol <- data() # reactive data must be read outside of a future
          content_type <- getQueryString()$content
          progress_bar <- Progress$new()
          progress_bar$set(message = "Building your shiny report", detail = "Patience, young grasshopper", value = 1)
          if (render_report_in_background) {
            result <- future({
              render_report_to_file(file, 
                                    input = gene_symbol, 
                                    type = content_type)
            })
            finally(result, function(){
              progress_bar$close()
            })
          } else {
            render_report_to_file(file,
                                  input = gene_symbol,
                                  type = content_type)
            progress_bar$close()
          }
        }
      )
    }
  )
} 
