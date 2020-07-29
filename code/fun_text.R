#this script contains functions for returning text, like summaries

# Gene Summary
summary_gene <- function(summary_table = gene_summary, gene_symbol, var = "approved_symbol") { #default so no error if empty
  gene_summary_var <- summary_table %>%
    dplyr::filter(approved_symbol == gene_symbol) %>% 
    dplyr::pull(var) #any column name
  return(gene_summary_var)
}

# Pathway Summary
summary_pathway <- function(summary_table = pathways, go_id, var = "pathway") {
  if (var == "data") {
    pathway_summary_var <- summary_table %>%
      dplyr::filter(go == go_id) %>% 
      unnest(data) %>% 
      dplyr::pull(gene) %>% 
      str_c(collapse = ", ")
  } else {
  pathway_summary_var <- summary_table %>%
    dplyr::filter(go == go_id) %>% 
    dplyr::pull(var)
  }
  return(pathway_summary_var)
}

# Gene list Summary
summary_gene_list <- function(summary_table = gene_summary, gene_list) {
  # Filter out invalid symbols for when a user edits "custom_gene_list" query parameter
  valid_gene_symbols <- summary_table %>%
    dplyr::filter(approved_symbol %in% gene_list) %>%
    dplyr::pull(approved_symbol) %>% 
    str_c(collapse = ", ")
  return(valid_gene_symbols)
}
