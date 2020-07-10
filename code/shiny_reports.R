downloadReportPanel <- function(id) {
  ns <- NS(id)
  tagList(
    h2("Report Generator"),    
    conditionalPanel(condition = paste0("input['", ns("submit"), "'] == 0"),
                     "Please enter your name and email address to download a report", 
                     textInput("first_name", "First Name", ""), 
                     textInput("last_name", "Last Name", ""), 
                     textInput("email", "Email Address", ""), 
                     actionButton(inputId = ns("submit"), 
                                  label = "Enter")),
    conditionalPanel(condition = paste0("input['", ns("submit"), "'] != 0"),
                     "To generate a report, click on the button below",
                     br(),
                     downloadButton(outputId = ns("report"), label = "Download report"))
  )
}

# Define the fields we want to save from the form
fields <- c("first_name", "last_name", "email")

save_data <- function(input) {
  # put variables in a data frame
  data <- data.frame(matrix(nrow=1,ncol=0))
  for (x in fields) {
    data[[x]] <- input[[x]]
  }
  data$submit_time <- date()
  
  # Create a unique file name
  file_name <- sprintf(
    "%s_%s.csv",
    as.integer(Sys.time()), 
    digest::digest(data) #gives it a unique name
  )
  
  #create dir
  directory_path <- here::here("user-data")
  dir.create(file.path(directory_path), showWarnings = FALSE)
  # Write the file to the local system as csv without column headers for ease of use
  write_csv(data, path=file.path(directory_path, file_name), col_names=FALSE)
}

downloadReportPanelServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {   
      observeEvent(input$submit, {
        save_data(input)
      })
      
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
                                    gene_symbol, 
                                    type = content_type,
                                    summary1 = gene_summary, 
                                    summary2 = pathways,
                                    cellbins_data = achilles, 
                                    expression_data = expression_join, 
                                    celldeps_data = achilles,
                                    mean = mean_virtual_achilles,
                                    cellanatogram_data = subcell,
                                    toptable_data = master_top_table, 
                                    bottomtable_data = master_bottom_table,
                                    enrichmenttop_data = master_positive, 
                                    enrichmentbottom_data = master_negative, 
                                    achilles_data = achilles)
            })
            finally(result, function(){
              progress_bar$close()
            })
          } else {
            render_report_to_file(file, 
                                  gene_symbol, 
                                  type = getQueryString()$content,
                                  summary1 = gene_summary, 
                                  summary2 = pathways,
                                  cellbins_data = achilles, 
                                  expression_data = expression_join, 
                                  celldeps_data = achilles,
                                  mean = mean_virtual_achilles,
                                  cellanatogram_data = subcell,
                                  toptable_data = master_top_table, 
                                  bottomtable_data = master_bottom_table,
                                  enrichmenttop_data = master_positive, 
                                  enrichmentbottom_data = master_negative, 
                                  achilles_data = achilles)
            progress_bar$close()
          }
        }
      )
    }
  )
} 
