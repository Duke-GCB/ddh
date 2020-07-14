#Cell Dep Plots -----
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
      output$text_cell_dep_plot <- renderText({paste0("Dependency plots generated for ", str_c(data(), collapse = ", "))})
      output$cell_deps <- renderPlotly({
        validate(
          need(data() %in% colnames(achilles), "No data found for this gene."))
        withProgress(message = 'Wait for it...', value = 1, {
          ggplotly(make_celldeps(achilles, expression_join, data(), mean_virtual_achilles), tooltip = "text")
        })
      })      
    }
  )
}

cellBinsPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(plotlyOutput(outputId = ns("cell_bins"))),
    tags$br(),
    fluidRow(tags$strong(plot_cellbins_title), plot_cellbins_legend)
  )
}

cellBinsPlotServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$cell_bins <- renderPlotly({
        validate(
          need(data() %in% colnames(achilles), "")) #""left blank
        ggplotly(make_cellbins(achilles, expression_join, data()), tooltip = c("text"))
      })
    }
  )
}

cellDepsLinPlot <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(plotlyOutput(outputId = ns("cell_deps_lin")),
    tags$br(),
    fluidRow(tags$strong(plot_celllin_title), plot_celllin_legend, 
             actionLink(inputId = ns("sublin_click"), "View sublineage plot below")), #add conditional panel for subtype
    tags$br(), 
    conditionalPanel(condition = paste0("input['", ns("sublin_click"), "'] != 0"), 
                     tags$br(),
                     fluidRow(plotlyOutput(outputId = ns("cell_deps_sublin")))))
  )
}

cellDepsLinPlotServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$cell_deps_lin <- renderPlotly({
        validate(
          need(data() %in% colnames(achilles), "No data found for this gene."))
        withProgress(message = 'Wait for it...', value = 1, {
          ggplotly(make_lineage(achilles, expression_join, data()), tooltip = "text")
        })
      })
      observeEvent(input$sublin_click, { #event to store the 'click'
      })
      output$cell_deps_sublin <- renderPlotly({
        validate(
          need(data() %in% colnames(achilles), "No data found for this gene."))
        withProgress(message = 'Wait for it...', value = 1, {
          ggplotly(make_sublineage(achilles, expression_join, data()), height = 1100, tooltip = "text")
        })
      })
    }
  )
}
