#Gene summary
geneSummaryText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("gene_summary_title"))),
    h4("Summary"),
    tags$dl(
      tags$dt("Gene"), tags$dd(textOutput(outputId = ns("gene_summary_approved_symbol"))),
      tags$dt("Name"), tags$dd(textOutput(outputId = ns("gene_summary_approved_name"))),
      tags$dt("aka"), tags$dd(textOutput(outputId = ns("gene_summary_aka"))),
      tags$dt("Entrez ID"), tags$dd(textOutput(outputId = ns("gene_summary_ncbi_gene_id"))),
      tags$dt("Gene Summary"), tags$dd(textOutput(outputId = ns("gene_summary_entrez_summary")))
    ),
  )
}

geneSummaryTextServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$gene_summary_title <- renderText({paste0(summary_gene(summary_table = gene_summary, gene_symbol = data(), var = "approved_symbol"), ": ", summary_gene(summary_table = gene_summary, gene_symbol = data(), var = "approved_name"))})
      output$gene_summary_approved_symbol <- renderText(summary_gene(summary_table = gene_summary, gene_symbol = data(), var = "approved_symbol"))
      output$gene_summary_approved_name <- renderText(summary_gene(summary_table = gene_summary, gene_symbol = data(), var = "approved_name"))
      output$gene_summary_aka <- renderText(summary_gene(summary_table = gene_summary, gene_symbol = data(), var = "aka"))
      output$gene_summary_ncbi_gene_id <- renderText(summary_gene(summary_table = gene_summary, gene_symbol = data(), var = "ncbi_gene_id"))
      output$gene_summary_entrez_summary <- renderText(summary_gene(summary_table = gene_summary, gene_symbol = data(), var = "entrez_summary"))
    })
}

#pathways
pathwaySummaryText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("pathway_summary_title"))),
    tags$dl(
      tags$dt("Genes"),
      tags$dd(textOutput(outputId = ns("pathway_summary_gene_symbols"))),
      tags$dt("Pathway Description"),
      tags$dd(textOutput(outputId = ns("pathway_summary_def")))
    ),
  )
}

pathwaySummaryTextServer <- function(id, pathway_go) {
  moduleServer(
    id,
    function(input, output, session) {
      output$pathway_summary_title <- renderText({paste0("Pathway: ", summary_pathway(summary_table = pathways, go_id = pathway_go, var = "pathway"), " (GO:", summary_pathway(summary_table = pathways, go_id = pathway_go, var = "go"), ")")})
      output$pathway_summary_gene_symbols <- renderText({summary_pathway(summary_table = pathways, go_id = pathway_go, var = "data")})
      output$pathway_summary_def <- renderText({summary_pathway(summary_table = pathways, go_id = pathway_go, var = "def")})
    }
  )
}

#gene_list
geneListSummaryText <- function (id) {
  ns <- NS(id)
  list(
    h3("Custom Gene List"),
    tags$dl(
      tags$dt("Genes"),
      tags$dd(textOutput(outputId = ns("custom_gene_list")))
    )
  )
}

geneListSummaryTextServer <- function(id, data) { #what is data here?
  moduleServer(
    id,
    function(input, output, session) {
      output$custom_gene_list <- renderText({summary_gene_list(summary_table = gene_summary, gene_list = data())})
    }
  )
}





