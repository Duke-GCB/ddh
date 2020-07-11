# QUERY SHINY FUNCTIONS ----
gene_query_result_row <- function(row) {
  gene_summary_row <- row$data
  title <- paste0(gene_summary_row["approved_symbol"], ": ", gene_summary_row["approved_name"])
  list(
    h4(
      tags$strong("Gene:"),
      tags$a(title, href=paste0("?show=gene&query_type=gene&symbol=", gene_summary_row["approved_symbol"]))
    ),
    div(tags$strong("Aka:"), gene_summary_row["aka"]),
    div(tags$strong("Entrez ID:"), gene_summary_row["ncbi_gene_id"]),
    hr()
  )
}

pathway_query_result_row <- function(row) {
  pathways_row <- row$data
  gene_symbols <- lapply(pathways_row$data, function(x) { paste(x$gene, collapse=', ') })
  title <- paste0(pathways_row$pathway, " (GO:", pathways_row$go, ")")
  list(
    h4(
      tags$strong("Pathway:"),
      tags$a(title, href=paste0("?show=pathway&query_type=pathway&go=", pathways_row$go))
    ),
    tags$dl(
      tags$dt("Genes"),
      tags$dd(gene_symbols),
    ),
    hr()
  )
}

gene_list_query_result_row <- function(row) {
  gene_summary_rows <- row$data
  title <- row$key
  
  known_gene_symbols <- gene_summary_rows %>% 
    filter(known == TRUE) %>%
    pull(approved_symbol)
  has_known_gene_symbols <- !is_empty(known_gene_symbols)
  
  unknown_gene_symbols <- gene_summary_rows %>% 
    filter(known == FALSE) %>%
    pull(approved_symbol)
  has_unknown_gene_symbols <- !is_empty(unknown_gene_symbols)
  
  known_gene_symbols_tags <- NULL
  if (has_known_gene_symbols) {
    gene_query_param <- paste0("custom_gene_list=", paste(known_gene_symbols, collapse=","))
    href <- paste0("?show=gene_list&query_type=custom_gene_list&", gene_query_param)
    known_gene_symbols_tags <- list(
      tags$h6("Known Gene Symbols"),
      tags$a(paste(known_gene_symbols, collapse=", "), href=href)
    )
  }
  
  unknown_gene_symbols_tags <- NULL
  if (has_unknown_gene_symbols) {
    unknown_gene_symbols_tags <- list(
      tags$h6("Unknown Gene Symbols"),
      tags$div(paste(unknown_gene_symbols, collapse=", "))
    )
  }
  
  list(
    h4(
      tags$strong("Custom Gene List")#,
      #tags$span(title)
    ),
    known_gene_symbols_tags,
    unknown_gene_symbols_tags,
    hr()
  )
}
