#Summary pages
#Gene summary
geneSummaryText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("gene_summary_title"))),
    h4("Summary"),
    tags$dl(
      tags$dt("Gene Summary"), tags$dd(textOutput(outputId = ns("gene_summary_entrez_summary")))
    ),
  )
}

geneSummaryTextServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$gene_summary_title <- renderText({paste0(summary_gene(summary_table = gene_summary, input = data(), var = "approved_symbol"), ": ", summary_gene(summary_table = gene_summary, input = data(), var = "approved_name"))})
      output$gene_summary_entrez_summary <- renderText(summary_gene(summary_table = gene_summary, input = data(), var = "entrez_summary"))
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

pathwaySummaryTextServer <- function(id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$pathway_summary_title <- renderText({paste0("Pathway: ", summary_pathway(summary_table = pathways, input = data(), var = "pathway"), " (GO:", summary_pathway(summary_table = pathways, input = data(), var = "go"), ")")})
      output$pathway_summary_gene_symbols <- renderText({summary_pathway(summary_table = pathways, input = data(), var = "data")})
      output$pathway_summary_def <- renderText({summary_pathway(summary_table = pathways, input = data(), var = "def")})
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
      output$custom_gene_list <- renderText({summary_gene_list(summary_table = gene_summary, input = data())})
    }
  )
}

#gene page
#this is the gene tab for a gene search
geneText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("gene_summary_title"))),
    h4("Summary"),
    tags$dl(
      tags$dt("Gene"), tags$dd(textOutput(outputId = ns("gene_summary_approved_symbol"))),
      tags$dt("Name"), tags$dd(textOutput(outputId = ns("gene_summary_approved_name"))),
      tags$dt("aka"), tags$dd(textOutput(outputId = ns("gene_summary_aka"))),
      tags$dt("Entrez ID"), tags$dd(htmlOutput(outputId = ns("ncbi_link"))),
      tags$dt("Gene Summary"), tags$dd(textOutput(outputId = ns("gene_summary_entrez_summary")))
    ),
  )
}

geneTextServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$gene_summary_title <- renderText({paste0(summary_gene(summary_table = gene_summary, input = data(), var = "approved_symbol"), ": ", summary_gene(summary_table = gene_summary, input = data(), var = "approved_name"))})
      output$gene_summary_approved_symbol <- renderText(summary_gene(summary_table = gene_summary, input = data(), var = "approved_symbol"))
      output$gene_summary_approved_name <- renderText(summary_gene(summary_table = gene_summary, input = data(), var = "approved_name"))
      output$gene_summary_aka <- renderText(summary_gene(summary_table = gene_summary, input = data(), var = "aka"))
      output$gene_summary_entrez_summary <- renderText(summary_gene(summary_table = gene_summary, input = data(), var = "entrez_summary"))
      output$ncbi_link <- renderText(paste0('<a href="https://www.ncbi.nlm.nih.gov/gene/?term=', summary_gene(summary_table = gene_summary, input = data(), var = "ncbi_gene_id"), '" target="_blank">', summary_gene(summary_table = gene_summary, input = data(), var = "ncbi_gene_id"),'</a>'))
      }
  )
}

#protein page
#this is the protein tab for a gene search
proteinText <- function (id) {
  ns <- NS(id)
  tagList(
    h3(textOutput(outputId = ns("protein_summary_title"))),
    h4("Summary"),
    tags$dl(
      tags$dt("Protein"), tags$dd(textOutput(outputId = ns("protein_summary_approved_symbol"))),
      tags$dt("Name"), tags$dd(textOutput(outputId = ns("protein_summary_approved_name"))),
      tags$dt("Protein Mass"), tags$dd(textOutput(outputId = ns("protein_summary_mass"))),
      tags$dt("Enzyme Commission", a(href = "https://en.wikipedia.org/wiki/Enzyme_Commission_number", img(src="link out_25.png", width="10", height="10"),  target="_blank")), tags$dd(htmlOutput(outputId = ns("ec_link"))),
      tags$dt("Uniprot ID"), tags$dd(htmlOutput(outputId = ns("uniprot_link"))),
      tags$dt("Protein Summary"), tags$dd(textOutput(outputId = ns("protein_summary_uniprot_summary"))), 
      tags$dt("Protein Sequence"), tags$dd(verbatimTextOutput(outputId = ns("protein_summary_seq"))),
    ),
  )
}

proteinTextServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
      output$protein_summary_title <- renderText({paste0(summary_protein(summary_table = proteins, input = data(), var = "gene_name"), ": ", summary_protein(summary_table = proteins, input = data(), var = "protein_name"))})
      output$protein_summary_approved_symbol <- renderText(summary_protein(summary_table = proteins, input = data(), var = "gene_name"))
      output$protein_summary_approved_name <- renderText(summary_protein(summary_table = proteins, input = data(), var = "protein_name"))
      output$ec_link <- renderText(
        if (is.na(summary_protein(summary_table = proteins, input = data(), var = "ec"))) {paste0('<a href="https://enzyme.expasy.org" target="_blank">', summary_protein(summary_table = proteins, input = data(), var = "ec"),'</a>')
        } else {paste0('<a href="https://enzyme.expasy.org/EC/', summary_protein(summary_table = proteins, input = data(), var = "ec"), '">', summary_protein(summary_table = proteins, input = data(), var = "ec"),'</a>')
        })
      output$protein_summary_uniprot_summary <- renderText(summary_protein(summary_table = proteins, input = data(), var = "function_cc"))
      output$uniprot_link <- renderText(paste0('<a href="https://www.uniprot.org/uniprot/', summary_protein(summary_table = proteins, input = data(), var = "uniprot_id"), '" target="_blank">', summary_protein(summary_table = proteins, input = data(), var = "uniprot_id"),'</a>'))
      output$protein_summary_mass <- renderText(paste0(summary_protein(summary_table = proteins, input = data(), var = "mass"), " kDa"))
      output$protein_summary_seq <- renderText(summary_protein(summary_table = proteins, input = data(), var = "sequence"))
    }
  )
}

#DUMMY
nameText <- function (id) {
    ns <- NS(id)
    "test"
}

nameTextServer <- function (id, data) {
  moduleServer(
    id,
    function(input, output, session) {
    }
  )
}

##TEMPLATE
# nameText <- function (id) {
#     ns <- NS(id)
# }
# 
# nameTextServer <- function (id, data) {
#   moduleServer(
#     id,
#     function(input, output, session) {
#     }
#   )
# }
