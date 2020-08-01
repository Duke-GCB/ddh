#this script contains functions for returning text, like summaries

# Gene Summary
summary_gene <- function(summary_table = gene_summary, input = gene_symbol, var = "approved_symbol") { #default so no error if empty
  gene_summary_var <- summary_table %>%
    dplyr::filter(approved_symbol == input) %>% 
    dplyr::pull(var) #any column name
  return(gene_summary_var)
}
# Pathway Summary
summary_pathway <- function(summary_table = pathways, input = go_id, var = "pathway") {
  if (var == "data") {
    pathway_summary_var <- summary_table %>%
      dplyr::filter(go == input) %>% 
      unnest(data) %>% 
      dplyr::pull(gene) %>% 
      str_c(collapse = ", ")
  } else {
  pathway_summary_var <- summary_table %>%
    dplyr::filter(go == input) %>% 
    dplyr::pull(var)
  }
  return(pathway_summary_var)
}

# Gene list Summary
summary_gene_list <- function(summary_table = gene_summary, input = gene_list) {
  # Filter out invalid symbols for when a user edits "custom_gene_list" query parameter
  valid_gene_symbols <- summary_table %>%
    dplyr::filter(approved_symbol %in% input) %>%
    dplyr::pull(approved_symbol) %>% 
    str_c(collapse = ", ")
  return(valid_gene_symbols)
}

#summary_gene(input = "SDHA", var = "entrez_summary")