downloadReportPanel <- function(id) {
  ns <- NS(id)
  tagList(
    h2("Report Generator"),
    textOutput(ns("help_message")),
    br(),
    conditionalPanel(condition = paste0("input['", ns("generate_report"), "'] == 0"),
      actionButton(ns("generate_report"), "Generate report")
    ),
    conditionalPanel(condition = paste0("output['", ns("report_zip_path"), "'] != ''"),
                     downloadButton(outputId = ns("report"), label = "Download report")
    )
  )
}

downloadReportPanelServer <- function(id, type, data) {
  moduleServer(
    id,
    function(input, output, session) {
      # reactive variable to hold path to the download zip
      report_zip_path <- reactiveVal("")
      # make download path an output so it can be used in a conditionalPanel
      output$report_zip_path <- renderText({
        report_zip_path()
      })
      # populate this variable even though it isn't displayed
      outputOptions(output, "report_zip_path", suspendWhenHidden = FALSE)

      # update help message based on status of report generation
      output$help_message <- renderText({
        if (input$generate_report == 0) {
          "To generate a report, click on the button below"
        } else {
          if (report_zip_path() != "") {
            "Report complete. Click on the button below to download."
          }
        }
      })

      # user clicks generate report save zip into temp_zip_dir
      observeEvent(input$generate_report, {
        # create a temporary directory and make it our working directory
        temp_zip_dir <- tempfile(pattern="tmpdir", tmpdir=here::here("report"))
        dir.create(temp_zip_dir)

        if(type == "gene"){
          filename <- paste0(data()$id, "_ddh.zip")
        } else if (type == "pathway") {
          filename <- paste0("go_", data()$id, "_ddh.zip")
        } else {
          filename <- paste0("custom_", paste(data()$gene_symbols, collapse="_"), "_ddh.zip")
        }
        filename <- file.path(temp_zip_dir, filename)

        data_values <- data() # reactive data must be read outside of a future
        progress_bar <- Progress$new()
        progress_bar$set(message = "Building your shiny report", detail = "Patience, young grasshopper", value = 1)
        if (render_report_in_background) {
          result <- future({
            render_report_to_file(data_values=data_values, file=filename)
            report_zip_path(filename)
          })
          finally(result, function(){
            progress_bar$close()
          })
        } else {
          render_report_to_file(data_values=data_values, file=filename)
          report_zip_path(filename)
          progress_bar$close()
        }
      })

      output$report <- downloadHandler(
        # create pdf report
        filename = function() {
          basename(report_zip_path())
        },
        content = function(file) {
          source_zip_path <- report_zip_path()
          message(file)
          message(source_zip_path)
          file.copy(source_zip_path, file)
        }
      )
    }
  )
} 
