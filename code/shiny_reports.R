downloadReportPanel <- function(id) {
  ns <- NS(id)
  tagList(
    h2("Report Generator"),    
                     "To generate a report, click on the button below",
                     br(),
                     downloadButton(outputId = ns("report"), label = "Download report")
  )
}

downloadReportPanelServer <- function(id, type, data) {
  moduleServer(
    id,
    function(input, output, session) {   
      output$report <- downloadHandler(
        # create pdf report
        filename = function() {
          if(type == "gene"){
            paste0(data()$id, "_ddh.zip")
          } else if (type == "pathway") {
            paste0("go_", data()$id, "_ddh.zip")
          } else {
            paste0("custom_", paste(data()$gene_symbols, collapse="_"), "_ddh.zip")
          }
        },
        content = function(file) {
          data_values <- data() # reactive data must be read outside of a future
          progress_bar <- Progress$new()
          progress_bar$set(message = "Building your shiny report", detail = "Patience, young grasshopper", value = 1)
          if (render_report_in_background) {
            result <- future({
              render_report_to_file(data_values=data_values, file=file)
            })
            finally(result, function(){
              progress_bar$close()
            })
          } else {
            render_report_to_file(data_values=data_values, file=file)
            progress_bar$close()
          }
        }
      )
    }
  )
} 
