# Plots -----
cellDependenciesPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_cell_dep_plot")))),
    fluidRow(plotlyOutput(outputId = ns("cell_deps"))),
    tags$br(),
    fluidRow(tags$strong(plot_celldeps_title), plot_celldeps_legend))
}

cellDependenciesPlotServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_cell_dep_plot <- renderText({paste0("Dependency plots generated for ", str_c(data()$gene_symbols, collapse = ", "))})
      output$cell_deps <- renderPlotly({
        validate(
          need(data()$gene_symbols %in% colnames(achilles), "No data found for this gene."))
        withProgress(message = 'Wait for it...', value = 1, {
          ggplotly(make_celldeps(achilles, expression_names, data()$gene_symbols, mean_virtual_achilles), tooltip = "text")
        })
      })      
    }
  )
}

cellBinsPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(plotOutput(outputId = ns("cell_bins"),  height = "auto")),
    tags$br(),
    fluidRow(tags$strong(plot_cellbins_title), plot_cellbins_legend)
  )
}

cellBinsPlotServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$cell_bins <- renderPlot({
        validate(
          need(data()$gene_symbols %in% colnames(achilles), "")) #""left blank
        make_cellbins(achilles, expression_names, data()$gene_symbols)
      },
      height = function() length(data()$gene_symbols) * 90 + 80)
    }
  )
}

cellDepsLinPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(plotOutput(outputId = ns("cell_deps_lin"),  height = "auto"),
    tags$br(),
    fluidRow(tags$strong(plot_celllin_title), plot_celllin_legend), 
    tags$br(),
    fluidRow(actionLink(inputId = ns("sublin_click"), "View sublineage plot below")),
    tags$br(), 
    conditionalPanel(condition = paste0("input['", ns("sublin_click"), "'] != 0"), 
                     tags$br(),
                     fluidRow(plotOutput(outputId = ns("cell_deps_sublin"),  height = "auto"))))
  )
}

cellDepsLinPlotServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$cell_deps_lin <- renderPlot({
        validate(
          need(data()$gene_symbols %in% colnames(achilles), "No data found for this gene."))
        withProgress(message = 'Wait for it...', value = 1, {
          make_lineage(achilles, expression_names, data()$gene_symbols)
        })
      },
      height = 550)
      observeEvent(input$sublin_click, { #event to store the 'click'
      })
      output$cell_deps_sublin <- renderPlot({
        validate(
          need(data()$gene_symbols %in% colnames(achilles), "No data found for this gene."))
        withProgress(message = 'Wait for it...', value = 1, {
          make_sublineage(achilles, expression_names, data()$gene_symbols)
        })
      },
      height = 1400)
    }
  )
}

cellAnatogramPlot <- function(id) {
  ns <- NS(id)
    fluidRow(plotOutput(outputId = ns("cellanatogram")))
}
  
cellAnatogramPlotServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$cellanatogram <- renderPlot({
        validate(
          need(data()$gene_symbols %in% subcell$gene_name, "No subcellular location data for this gene."))
        make_cellanatogram(subcell, data()$gene_symbols)
      })
    }
  )
}

cellExpressionPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(h4(textOutput(ns("text_cell_exp_plot")))),
    fluidRow(plotlyOutput(outputId = ns("cell_exp"))),
  )
}

cellExpressionPlotServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$text_cell_exp_plot <- renderText({paste0("Expression plots generated for ", str_c(data()$gene_symbols, collapse = ", "))})
      output$cell_exp <- renderPlotly({
        validate(
          need(data()$gene_symbols %in% colnames(expression), "")) #""left blank
        ggplotly(make_cellexpression(gene_symbol = data()$gene_symbols), tooltip = c("text"))
      })
    }
  )
}
