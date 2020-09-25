library(tidyverse)

make_top_table <- function(toptable_data = master_top_table, gene_symbol) {
  toptable_data %>%
    dplyr::filter(fav_gene %in% gene_symbol) %>%
    tidyr::unnest(data) %>%
    dplyr::arrange(desc(r2)) %>%
    dplyr::rename("Query" = "fav_gene", "Gene" = "gene", "Name" = "name", "R^2" = "r2", "Z Score" = "z_score", "Co-publication Count" = "concept_count", "Co-publication Index" = "concept_index")
}

censor <- function(top_table, censor_data = censor_genes, choice = FALSE, greater_than) {
  if(choice == TRUE){
    censor_data <- censor_data %>%
      dplyr::filter(num_sim > greater_than)
    
    censored_table <- top_table %>%
      dplyr::filter(!Gene %in% censor_data$genes)
    
    return(censored_table)
  }
  return(top_table)
}

make_bottom_table <- function(bottomtable_data = master_bottom_table, gene_symbol) {
  bottomtable_data %>%
    dplyr::filter(fav_gene %in% gene_symbol) %>%
    tidyr::unnest(data) %>%
    dplyr::arrange(r2) %>%
    dplyr::rename("Query" = "fav_gene", "Gene" = "gene", "Name" = "name", "R^2" = "r2", "Z Score" = "z_score", "Co-publication Count" = "concept_count", "Co-publication Index" = "concept_index")
}

make_enrichment_top <- function(enrichmenttop_data, gene_symbol) { #master_positive
  enrichmenttop_data %>%
    dplyr::filter(fav_gene %in% gene_symbol) %>%
    tidyr::unnest(data) %>%
    dplyr::select(fav_gene, enrichr, Term, Overlap, Adjusted.P.value, Combined.Score, Genes) %>%
    dplyr::arrange(Adjusted.P.value) %>%
    dplyr::rename("Query" = "fav_gene", "Gene Set" = "enrichr", "Gene List" = "Term", "Adjusted p-value" = "Adjusted.P.value", "Combined Score" = "Combined.Score") #"Overlap", "Genes"
}

make_enrichment_bottom <- function(enrichmentbottom_data, gene_symbol) { #master_negative
  enrichmentbottom_data %>%
    dplyr::filter(fav_gene %in% gene_symbol) %>%
    tidyr::unnest(data) %>%
    dplyr::select(fav_gene, enrichr, Term, Overlap, Adjusted.P.value, Combined.Score, Genes) %>%
    dplyr::arrange(Adjusted.P.value) %>%
    dplyr::rename("Query" = "fav_gene", "Gene Set" = "enrichr", "Gene List" = "Term", "Adjusted p-value" = "Adjusted.P.value", "Combined Score" = "Combined.Score") #"Overlap", "Genes"
}

make_achilles_table <- function(achilles_data, expression_data, gene_symbol) { #achilles, expression_join, gene_symbol
  target_achilles <- achilles_data %>%
    dplyr::select(X1, any_of(gene_symbol)) %>%
    dplyr::left_join(expression_data, by = "X1") %>%
    dplyr::select(-X1) %>%
    dplyr::select(cell_line, lineage, lineage_subtype, everything()) %>% 
    dplyr::mutate_if(is.numeric, ~round(., digits = 3)) %>% 
    dplyr::rename("Cell Line" = "cell_line", "Lineage" = "lineage", "Subtype" = "lineage_subtype") %>% 
    dplyr::arrange(.[[4]])
  return(target_achilles)
}

make_pathway_table <- function(table_name = pathways) { #this makes a summary table of pathways for browsing
  pathway_table <- table_name %>% 
    mutate(genes = map_chr(table_name$data, function(.x) {
      y <- unlist(.x, use.names = FALSE)
      str_c(y, collapse = ', ')
    })) %>% 
    select(pathway, go, genes)
  return(pathway_table)
}

word_starts_with_regex <- function(query_str) {
  regex(paste0('\\b', query_str), ignore_case=TRUE)
}

search_tables <- function(gene_summary, pathways, expression_names, query_str) {
  if (grepl(',', query_str)) {
    custom_list_search_tables(gene_summary, expression_names, query_str)
  } else {
    regular_search_tables(gene_summary, pathways, expression_names, query_str)
  }
}

custom_list_search_tables <- function(gene_summary, expression_names, query_str) {
  bind_rows(
    custom_gene_list_search_tables(gene_summary, query_str),
    custom_cell_line_list_search_tables(expression_names, query_str)
  )
}

regular_search_tables <- function(gene_summary, pathways, expression_names, query_str, limit_rows=10) {
  gene_data_search_result <- search_gene_data(gene_summary, pathways, query_str, limit_rows)
  cell_line_search_result <- search_cell_line_data(expression_names, query_str, limit_rows)
  bind_rows(
    gene_data_search_result,
    cell_line_search_result
  )
}

search_cell_line_data <- function(expression_names, query_str, limit_rows) {
  word_starts_with_query_str <- word_starts_with_regex(query_str)

  # find cell lines that start with query_str
  cell_line_rows <- expression_names %>%
    filter(str_detect(cell_line, word_starts_with_query_str)) %>%
    mutate(length = str_count(.[[1]])) %>%
    arrange(length) %>%
    head(limit_rows) %>%
    select(-length)

  # group cell lines into generic grouped format
  cell_line_data <- cell_line_rows %>%
    head(limit_rows) %>%
    mutate(key = cell_line) %>%
    add_column(query_type='cell') %>%
    group_by(key, query_type) %>%
    nest()

  # find lineage that start with query_str
  lineage_rows <- expression_names %>%
    filter(str_detect(lineage, word_starts_with_query_str)) %>%
    group_by(lineage) %>%
    group_nest() %>%
    mutate(length = str_count(.[[1]])) %>%
    arrange(length) %>%
    head(limit_rows) %>%
    select(-length)

  # group lineage data into generic grouped format
  lineage_data <- lineage_rows %>%
    head(limit_rows) %>%
    mutate(key = lineage,
           title = lineage) %>%
    add_column(query_type='lineage') %>%
    group_by(key, title, query_type) %>%
    nest()

  # find lineage subtype that start with query_str
  lineage_subtype_rows <- expression_names %>%
    filter(str_detect(lineage_subtype, word_starts_with_query_str)) %>%
    group_by(lineage_subtype) %>%
    group_nest() %>%
    mutate(length = str_count(.[[1]])) %>%
    arrange(length) %>%
    head(limit_rows) %>%
    select(-length)

  # group lineage data into generic grouped format
  lineage_subtype_data <- lineage_subtype_rows %>%
    head(limit_rows) %>%
    mutate(key = lineage_subtype,
           title = lineage_subtype) %>%
    add_column(query_type='lineage_subtype') %>%
    group_by(key, title, query_type) %>%
    nest()

  bind_rows(cell_line_data, lineage_data, lineage_subtype_data)
}

search_gene_data <- function(gene_summary, pathways, query_str, limit_rows) {
  word_starts_with_query_str <- word_starts_with_regex(query_str)

  # find pathway data
  pathways_data_title <- pathways %>%
    filter(str_detect(pathway, word_starts_with_query_str)) %>%
    mutate(length = str_count(.[[1]])) %>% 
    arrange(length) %>% 
    head(limit_rows) %>%
    select(-length)

  # find go data
  pathways_data_go <- pathways %>%
    filter(str_detect(go, word_starts_with_query_str)) %>%
    head(limit_rows)

  # nest pathways data underneath generic key, title, and query_type columns
  pathways_data <- unique(bind_rows(pathways_data_title, pathways_data_go)) %>%
    head(limit_rows) %>%
    mutate(key = go,
           title = pathway) %>%
    add_column(query_type='pathway') %>%
    group_by(key, title, query_type) %>%
    nest()

  # find genes most specific
  genes_data_symbol <- gene_summary %>%
    filter(str_detect(approved_symbol, word_starts_with_query_str)) %>%
    mutate(length = str_count(.[[1]])) %>% 
    arrange(length) %>% 
    head(limit_rows) %>%
    select(-length)

  # find genes most likely alternative
  genes_data_aka <- gene_summary %>%
    filter(str_detect(aka, word_starts_with_query_str)) %>%
    mutate(length = str_count(.[[1]])) %>% 
    arrange(length) %>% 
    head(limit_rows) %>%
    select(-length)

  # find genes most generic
  genes_data_name <- gene_summary %>%
    filter(str_detect(approved_name, word_starts_with_query_str)) %>%
    mutate(length = str_count(.[[1]])) %>% 
    arrange(length) %>% 
    head(limit_rows) %>%
    select(-length)

  # nest gene data underneath generic key, title, and query_type columns
  genes_data <- unique(bind_rows(genes_data_symbol, genes_data_aka, genes_data_name)) %>%
    head(limit_rows) %>%
    mutate(key = approved_symbol) %>%
    add_column(query_type='gene') %>%
    group_by(key, query_type) %>%
    nest()

  bind_rows(genes_data, pathways_data)
}

custom_gene_list_search_tables <- function(gene_summary, query_str) {
  # Search in gene_summary for matching approved_symbols
  # Returns df with columns:
  # - key - query_str passed in
  # - query_type - 'gene_list'
  # - data - gene_summary row for the specified key each gene symbol, when not found only the approved_symbol will be populated
  # create a df containing valid gene summary rows and just the approved_symbol filled in for unknown gene symbols
  query_gene_symbols <- c(str_split(query_str, "\\s*,\\s*", simplify = TRUE))
  gene_symbol_with_known <- query_gene_symbols %>%
    map_dfr(query_symbol_in_gene_summary, gene_summary=gene_summary)
  if(any(gene_symbol_with_known$known)) {
    gene_symbol_with_known %>%
      add_column(key=query_str) %>%
      add_column(query_type='gene_list') %>%
      group_by(key, query_type) %>%
      nest()
  } else {
    tibble()
  }
}

query_symbol_in_gene_summary <- function(gene_symbol, gene_summary) {
  # Searches for exact match on approved_symbol in gene_summary, when not found returns row with only approved_symbol filled in
  matches_gene_symbol_ignore_case <- regex(paste0('^', gene_symbol, '$'), ignore_case = TRUE)
  gene_summary_row <- gene_summary %>%
    filter(str_detect(approved_symbol, matches_gene_symbol_ignore_case))
  if (nrow(gene_summary_row) > 0) {
    gene_summary_row  %>%
      add_column(known=TRUE)
  } else {
    tibble(approved_symbol=gene_symbol, known=FALSE)
  }
}

custom_cell_line_list_search_tables <- function(expression_names, query_str) {
  cell_line_symbols <- c(str_split(query_str, "\\s*,\\s*", simplify = TRUE))
  cell_line_symbols_with_known <- cell_line_symbols %>%
    map_dfr(query_cell_line_in_expression_names, expression_names=expression_names)
  if(any(cell_line_symbols_with_known$known)) {
    cell_line_symbols_with_known %>%
      add_column(key=query_str) %>%
      add_column(query_type='cell_list') %>%
      group_by(key, query_type) %>%
      nest()
  } else {
    tibble()
  }
}

query_cell_line_in_expression_names <- function(cell_line, expression_names) {
  matches_cell_line_ignore_case <- regex(paste0('^', cell_line, '$'), ignore_case = TRUE)
  expression_names_row <- expression_names %>%
    filter(str_detect(cell_line, matches_cell_line_ignore_case))
  if (nrow(expression_names_row) > 0) {
    expression_names_row  %>%
      add_column(known=TRUE)
  } else {
    tibble(cell_line=cell_line, known=FALSE)
  }
}

make_cellanatogram_table <- function(cellanatogram_data = subcell, gene_symbol) {
  cellanatogram_data %>% 
    filter_all(any_vars(gene_name %in% gene_symbol)) %>% 
    filter(!is.na(type)) %>%
    add_count(main_location) %>% 
    transmute(Gene = gene_name, Reliability = reliability, Location = main_location, Count = as_factor(n)) %>% 
    arrange(desc(Count))
}

make_expression_table <- function(expression_data = expression, expression_join = expression_names, gene_symbol) { #you are so slow
  expression_data %>% 
    dplyr::select(is.character, any_of(gene_symbol)) %>% 
    dplyr::left_join(expression_join, by = "X1") %>%
    dplyr::select(-X1) %>%
    dplyr::select(cell_line, lineage, lineage_subtype, everything()) %>% 
    dplyr::mutate_if(is.numeric, ~round(., digits = 3)) %>% 
    dplyr::rename("Cell Line" = "cell_line", "Lineage" = "lineage", "Subtype" = "lineage_subtype") %>% 
    dplyr::arrange(desc(.[[4]]))
}
