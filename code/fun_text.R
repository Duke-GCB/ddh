#this script contains functions for returning text, like summaries

# Gene Summary
summary_gene <- function(summary_table = gene_summary, input = list(), var = "approved_symbol") { #default so no error if empty, but this pulls the var out of the df
  gene_summary_var <- summary_table %>%
    dplyr::filter(approved_symbol == input$id) %>%
    dplyr::pull(var) #any column name
  return(gene_summary_var)
}
# Pathway Summary
summary_pathway <- function(summary_table = pathways, input = list(), var = "pathway") {
  if (var == "data") {
    pathway_summary_var <- summary_table %>%
      dplyr::filter(go == input$id) %>%
      unnest(data) %>% 
      dplyr::pull(gene) %>% 
      str_c(collapse = ", ")
  } else {
  pathway_summary_var <- summary_table %>%
    dplyr::filter(go == input$id) %>%
    dplyr::pull(var)
  }
  return(pathway_summary_var)
}

# Gene list Summary
summary_gene_list <- function(summary_table = gene_summary, input = list()) {
  # Filter out invalid symbols for when a user edits "custom_gene_list" query parameter
  valid_gene_symbols <- summary_table %>%
    dplyr::filter(approved_symbol %in% input$gene_symbols) %>%
    dplyr::pull(approved_symbol) %>% 
    str_c(collapse = ", ")
  return(valid_gene_symbols)
}

#summary_gene(input = "SDHA", var = "entrez_summary")

# protein summary
#this is a master function to pull data out of the proteins df
summary_protein <- function(summary_table = proteins, input = list(), var = "gene_name") { #default so no error if empty, but this pulls the var out of the df
  protein_summary_var <- summary_table %>%
    dplyr::filter(gene_name %in% input$id) %>%
    dplyr::pull(var) #any column name
  return(protein_summary_var)
}
#summary_protein(input = list(id = "GSS"), var = "gene_name")
#summary_protein(input = list(id = "GSS"), var = "protein_name")

