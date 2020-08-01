library(tidyverse)

make_summary <- function(input, type) { #do I need to carry over summary table vars?
  if(type == "gene"){
    summary_table <- tibble(
      identifier = summary_gene(summary_table = gene_summary, input, var = "approved_symbol"),
      name = summary_gene(summary_table = gene_summary, input, var = "approved_name"),
      summary = summary_gene(summary_table = gene_summary, input, var = "entrez_summary"))
  } else if (type == "pathway") {
    summary_table <- tibble(
      identifier = paste0("GO:", summary_pathway(summary_table = pathways, input, var = "go")),
      name = summary_pathway(summary_table = pathways, input, var = "pathway"),
      summary = summary_pathway(summary_table = pathways, input, var = "def"))
  } else { #gene_list
    summary_table <- tibble(
      identifier = c("custom gene list"),
      name = summary_gene_list(summary_table = gene_summary, input),
      summary = c("User defined gene input list"))
  }
  return(summary_table)
}

#tests
#make_summary(input = "GSS", type = "gene")
#make_summary(input = "0060148", type = "pathway")
#make_summary(input = c("SDHA", "SDHB"), type = "gene_list")

#render in temp dir replaces usual render function
render_rmarkdown_in_tempdir <- function(rmd_path, output_file, envir = parent.frame()) {
  # The rmd_path variable must be an absolute path.
  
  # make sure the base report directory exists
  report_base_dir = here::here("report")
  if (!file.exists(report_base_dir)) {
    dir.create(report_base_dir)
  }
  # determine the filename of the Rmd file we will use for rendering
  rmd_filename <- basename(rmd_path)
  # create a temporary directory and make it our working directory
  temp_dir <- tempfile(pattern="tmpdir", tmpdir=report_base_dir)
  dir.create(temp_dir)
  owd <- setwd(temp_dir)
  on.exit(setwd(owd))
  on.exit(unlink(temp_dir, recursive = TRUE))
  # copy the Rmd file into our temporary(current) directory
  file.copy(rmd_path, rmd_filename, overwrite = TRUE)
  rmarkdown::render(rmd_filename, output_file = output_file, envir = envir)
}

#specific instructions to render reports based on query type and report template
render_gene_report <- function(input, type, output_file) {
  if(type == "gene"){
    gene_symbol <- input
  } else if (type == "pathway") {
    gene_symbol <- pathways %>% 
      filter(go %in% input) %>% 
      unnest(data) %>% 
      pull(gene)
  } else if (type == "gene_list") {
    gene_symbol <- input
  } else {
    stop("delcare your type!")
  }
  num <- length(achilles$X1)
  summary <- make_summary(input, type)
  cellanatogram <- make_cellanatogram(cellanatogram_data = subcell, gene_symbol)
  cellanatogram_table <- make_cellanatogram_table(cellanatogram_data = subcell, gene_symbol)
  p1 <- make_celldeps(celldeps_data = achilles, expression_data = expression_join, gene_symbol, mean = mean_virtual_achilles)
  p2 <- make_cellbins(cellbins_data = achilles, expression_data = expression_join, gene_symbol)
  p3 <- make_lineage(celldeps_data = achilles, expression_data = expression_join, gene_symbol)
  p4 <- make_sublineage(celldeps_data = achilles, expression_data = expression_join, gene_symbol)
  target_achilles_bottom <- make_achilles_table(achilles_data = achilles, expression_data = expression_join, gene_symbol) %>% head(10)
  target_achilles_top <- make_achilles_table(achilles_data = achilles, expression_data = expression_join, gene_symbol) %>% tail(10)
  dep_top <- make_top_table(toptable_data = master_top_table, gene_symbol)
  flat_top_complete <- make_enrichment_top(enrichmenttop_data = master_positive, gene_symbol)
  dep_bottom <- make_bottom_table(bottomtable_data = master_bottom_table, gene_symbol)
  flat_bottom_complete <- make_enrichment_bottom(enrichmentbottom_data = master_negative, gene_symbol)
  graph_report <- make_graph_report(toptable_data = master_top_table, bottomtable_data = master_bottom_table, gene_symbol)
  render_rmarkdown_in_tempdir(here::here("code", "report_gene.Rmd"), output_file)
}

#render_gene_report(input = "SST", type = "gene", output_file = "sst.pdf")
#render_gene_report(input = "0060148", type = "pathway")
#render_gene_report(input = c("GSS", "SST"), type = "gene_list")


#logic to matching query type to rendered content
render_report_to_file <- function(input,
                                  type, 
                                  file) {
  if (type == "gene") {
    render_gene_report(input, type, output_file = file)
  } else {
    stop("no report for you!")
  }
}

#render_report_to_file(input = "GSS", type = "gene", file = "gss_trial.pdf")

