#Cell Dep Table -----
cellDependenciesTable <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_cell_dep_table")))),
    fluidRow(dataTableOutput(outputId = ns("target_achilles"))))
}

cellDependenciesTableServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_cell_dep_table <- renderText({paste0("Dependency table generated for ", str_c(data(), collapse = ", "))})
      output$target_achilles <- DT::renderDataTable({
        validate(
          need(data() %in% colnames(achilles), "No data found for this gene."))
        make_achilles_table(achilles, expression_join, data())
      })
    }
  )
}


similarGenesTable <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_dep_top"))),
    fluidRow(
      column(6, checkboxGroupInput(inputId = ns("vars_dep_top"), 
                                   "Select columns:",
                                   c("R^2", "Z Score", "Co-publication Count", "Co-publication Index"), 
                                   selected = c("Z Score", "Co-publication Count"), 
                                   inline = TRUE)),
      column(6, fluidRow(sliderInput(inputId = ns("num_sim_genes"),
                                     "Censor genes with more than n associations:",
                                     min = 100,
                                     max = 1000,
                                     value = 1000,
                                     step = 100)), 
             fluidRow(column(3, actionButton(inputId = ns("censor"), "Submit")),
                      column(3, actionButton(inputId = ns("reset"), "Reset"))))
    ),
    hr(),
    fluidRow(dataTableOutput(outputId = ns("dep_top"))))
  )
}

similarGenesTableServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) { 
      output$text_dep_top <- renderText({paste0("Genes with similar dependencies as ", str_c(data(), collapse = ", "))})      
      output$dep_top <- DT::renderDataTable({
        validate(
          need(data() %in% master_top_table$fav_gene, "No data found for this gene."))
        DT::datatable(
          make_top_table(master_top_table, data()) %>% 
            dplyr::mutate(link = paste0("<center><a href='?show=detail&content=gene&symbol=", Gene,"'>", img(src="link out_25.png", width="10", height="10"),"</a></center>")) %>% 
            dplyr::select("Query", "Gene", "Gene \nLink" = "link", "Name", input$vars_dep_top) %>%
            censor(censor_genes, censor_status$choice, censor_status$num_sim_genes),
          escape = FALSE,
          options = list(pageLength = 25))
      })
      #censor reactive values
      censor_status <- reactiveValues(choice = FALSE, 
                                      num_sim_genes = 1000)
      observeEvent(input$censor, {
        censor_status$choice <- TRUE
        censor_status$num_sim_genes <- input$num_sim_genes
      })
      
      observeEvent(input$reset, {
        censor_status$choice <- FALSE
        censor_status$num_sim_genes <- 1000
        updateSliderInput(session, inputId = "num_sim_genes", value = 1000)
      })      
    }
  )
}


similarPathwaysTable <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_pos_enrich")))),
    fluidRow(dataTableOutput(outputId = ns("pos_enrich")))
  )
}

similarPathwaysTableServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) { 
      output$text_pos_enrich <- renderText({paste0("Pathways of genes with similar dependencies as ", str_c(data(), collapse = ", "))})
      output$pos_enrich <- DT::renderDataTable({
        validate(
          need(data() %in% master_positive$fav_gene, "No data found for this gene."))
        DT::datatable(
          make_enrichment_top(master_positive, data()),
          options = list(pageLength = 25))
      })      
    }
  )
}


dissimilarGenesTable <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_dep_bottom")))),
    fluidRow(checkboxGroupInput(inputId = ns("vars_dep_bottom"), 
                                "Select columns:",
                                c("R^2", "Z Score", "Co-publication Count", "Co-publication Index"), 
                                selected = c("Z Score", "Co-publication Count"), 
                                inline = TRUE)),
    fluidRow(dataTableOutput(outputId = ns("dep_bottom")))
  )
}

dissimilarGenesTableServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) { 
      output$text_dep_bottom <- renderText({paste0("Genes with inverse dependencies as ", str_c(data(), collapse = ", "))})      
      output$dep_bottom <- DT::renderDataTable({
        validate(
          need(data() %in% master_bottom_table$fav_gene, "No data found for this gene."))
        DT::datatable(
          make_bottom_table(master_bottom_table, data()) %>% 
            dplyr::mutate(link = paste0("<center><a href='?show=detail&content=gene&symbol=", Gene,"'>", img(src="link out_25.png", width="10", height="10"),"</a></center>")) %>% 
            dplyr::select("Query", "Gene", "Gene \nLink" = "link", "Name", input$vars_dep_bottom),
          escape = FALSE,
          options = list(pageLength = 25))
      })      
    }
  )
}


dissimilarPathwaysTable <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_neg_enrich")))),
    fluidRow(dataTableOutput(outputId = ns("neg_enrich")))
  )
}

dissimilarPathwaysTableServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) { 
      output$text_neg_enrich <- renderText({paste0("Pathways of genes with inverse dependencies as ", str_c(data(), collapse = ", "))})
      output$neg_enrich <- DT::renderDataTable({
        validate(
          need(data() %in% master_negative$fav_gene, "No data found for this gene."))
        DT::datatable(
          make_enrichment_bottom(master_negative, data()),
          options = list(pageLength = 25))
      })      
    }
  )
}


