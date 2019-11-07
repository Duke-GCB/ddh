library(tidyverse)
library(shiny)
library(RPostgreSQL)
library(feather)
library(shiny)
library(here)

drv <- dbDriver("PostgreSQL")
achilles_path <- "../data/19Q3_achilles.feather"
achilles <- read_feather(achilles_path)
achilles_long <- achilles %>% gather("gene", "dep_score", -X1)
expression_id <- read_feather("../data/19Q3_expression_id.feather")

expression_join <- expression_id %>% 
  rename(X1 = dep_map_id) %>% 
  select(X1, stripped_cell_line_name, lineage)

connectToDb <- function() {
    dbConnect(drv, user='postgres', dbname='depmap', host='localhost')
}

query_gene_summary <- function(con, gene_symbol) {
    sql <- "SELECT *
            FROM gene_summary
            WHERE approved_symbol = $1"
    gene_summary_rows <- dbGetQuery(con, sql, gene_symbol)
    if (nrow(gene_summary_rows) > 0) {
        gene_summary_row <- gene_summary_rows[1,]
    } else {
        NULL
    }
}

geneSummaryTagList <- function(gene_summary) {
    title <- paste0(gene_summary$approved_symbol, ": ", gene_summary$approved_name)
    tagList(
        h3(title),
        h4("Summary"),
        tags$dl(
            tags$dt("Gene"), tags$dd(gene_summary$approved_symbol),
            tags$dt("Name"), tags$dd(gene_summary$approved_name),
            tags$dt("aka"), tags$dd(gene_summary$aka),
            tags$dt("Entrez ID"), tags$dd(gene_summary$ncbi_gene_id),
        ),
        downloadButton("report", "Generate report"),
    )
}

render_report <- function(file, fav_gene, fav_gene_summary) {
    target_achilles <- achilles_long %>% 
      filter(gene == fav_gene) %>% 
      left_join(expression_join, by = "X1") %>% 
      select(stripped_cell_line_name, lineage, dep_score)
    dep_plot1 <- ggplot(target_achilles) +
      geom_histogram(aes(x = dep_score), binwidth = 0.25, color = "lightgray") +
      labs(x = "Dependency Score (binned)") + 
      theme_light()    
    
    # Copy the report file to a temporary directory before processing it, in
    # case we don't have write permissions to the current working dir (which
    # can happen when deployed).
    tempReport <- file.path(tempdir(), "report.Rmd")
    file.copy("../code/web_report.Rmd", tempReport, overwrite = TRUE)
    
    num_cell_lines <- length(achilles$X1)
    fav_gene <- fav_gene
    
    rmarkdown::render(tempReport, output_file = file)
}

ui <- fluidPage(
    titlePanel("Depmap"),
    sidebarLayout(
        sidebarPanel(
            textInput("gene_symbol", "Enter gene symbol", "", placeholder='BRCA1')
        ),
        mainPanel(
            uiOutput("ui")
        )
    )
)

server <- function(input, output, session) {
    observe({
        query <- parseQueryString(session$clientData$url_search)
        query_gene <- query[['gene']]
        if (!is.null(query_gene)) {
            updateTextInput(session, "gene_symbol", value=query_gene)
        }
    })
    output$ui <- renderUI({
        result <- tagList()
        if (input$gene_symbol != '') {
            con <- connectToDb()
            gene_summary <- query_gene_summary(con, input$gene_symbol)
            dbDisconnect(con)            
            if (is.null(gene_summary)) {
                result <- tagList(
                    h4(paste0("Gene symbol not found ", input$gene_symbol))
                )
            } else {
                title <- paste0(gene_summary$approved_symbol, ": ", gene_summary$approved_name)
                result <- geneSummaryTagList(gene_summary)
            }
        }
        result
    })
    output$report <- downloadHandler(
        filename = paste0(input$gene_symbol, ".pdf"),
        content = function(file) {
            con <- connectToDb()
            gene_summary <- query_gene_summary(con, input$gene_symbol)
            dbDisconnect(con)
            if (!is.null(gene_summary)) {
                render_report(file, input$gene_symbol, gene_summary)
            }
        }
    )
}

shinyApp(ui, server)
