library(tidyverse)

make_summary <- function(data_values) { #do I need to carry over summary table vars?
  if(data_values$type == "gene"){
    summary_table <- tibble(
      identifier = summary_gene(summary_table = gene_summary, data_values, var = "approved_symbol"),
      name = summary_gene(summary_table = gene_summary, data_values, var = "approved_name"),
      summary = summary_gene(summary_table = gene_summary, data_values, var = "entrez_summary"))
  } else if (data_values$type == "pathway") {
    summary_table <- tibble(
      identifier = paste0("GO:", summary_pathway(summary_table = pathways, data_values, var = "go")),
      name = summary_pathway(summary_table = pathways, data_values, var = "pathway"),
      summary = summary_pathway(summary_table = pathways, data_values, var = "def"))
  } else { #gene_list
    summary_table <- tibble(
      identifier = c("custom gene list"),
      name = summary_gene_list(summary_table = gene_summary, data_values),
      summary = c("User defined gene input list"))
  }
  return(summary_table)
}

#tests
#make_summary(data_values = list(id="GSS", type="gene", gene_symbols=c("GSS"))
#make_summary(data_values = list(id="0060148", type="pathway", gene_symbols=c("SDHA", "SDHB"))
#make_summary(data_values = list(id="SDHA,SDHB", gene_symbols=c("SDHA", "SDHB"), type="gene_list")

#render in temp dir replaces usual render function
render_rmarkdown_in_tempdir <- function(data_values, rmd_path, output_file, envir = parent.frame()) {
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

  #good file names
  good_file_name <- data_values$id
  if (data_values$type == "gene_list") {
    good_file_name <- paste0("custom_", paste(data_values$gene_symbols, collapse="_"))
  }
  #zip
  output_pdf_filename <- paste0(good_file_name, "_report.pdf")
  zip_filenames <- c(output_pdf_filename)
  rmarkdown::render(rmd_filename, output_file = output_pdf_filename, envir = envir)
  # get the names of all the items included for rendering
  for (name in names(envir)) {
    env_item = envir[[name]]
    # if the env_item is a plot
    if ("plot_env" %in% names(env_item)) {
      # save the plot to a png
      plot_filename <- paste0(good_file_name, "_", name, ".png")
      ggsave(plot_filename, env_item)
      # include the plot png in the zip download
      zip_filenames <- append(zip_filenames, plot_filename)
    }
  }  
  zip(zipfile = output_file, files = zip_filenames)
}

#specific instructions to render reports based on query type and report template
render_gene_report <- function(data_values, output_file) {
  if(data_values$type == "gene" | data_values$type == "pathway" | data_values$type == "gene_list") {
    gene_symbol <- data_values$gene_symbols
  } else {
    stop("delcare your type!")
  }
  num <- length(achilles$X1)
  summary <- make_summary(data_values)
  cellanatogram <- make_cellanatogram(cellanatogram_data = subcell, gene_symbol)
  cellanatogram_table <- make_cellanatogram_table(cellanatogram_data = subcell, gene_symbol)
  celldeps <- make_celldeps(celldeps_data = achilles, expression_data = expression_join, gene_symbol, mean = mean_virtual_achilles)
  cellbins <- make_cellbins(cellbins_data = achilles, expression_data = expression_join, gene_symbol)
  lineage <- make_lineage(celldeps_data = achilles, expression_data = expression_join, gene_symbol)
  sublineage <- make_sublineage(celldeps_data = achilles, expression_data = expression_join, gene_symbol)
  target_achilles_bottom <- make_achilles_table(achilles_data = achilles, expression_data = expression_join, gene_symbol) %>% head(10)
  target_achilles_top <- make_achilles_table(achilles_data = achilles, expression_data = expression_join, gene_symbol) %>% tail(10)
  dep_top <- make_top_table(toptable_data = master_top_table, gene_symbol)
  flat_top_complete <- make_enrichment_top(enrichmenttop_data = master_positive, gene_symbol)
  dep_bottom <- make_bottom_table(bottomtable_data = master_bottom_table, gene_symbol)
  flat_bottom_complete <- make_enrichment_bottom(enrichmentbottom_data = master_negative, gene_symbol)
  graph <- make_graph_report(toptable_data = master_top_table, bottomtable_data = master_bottom_table, gene_symbol)
  render_rmarkdown_in_tempdir(data_values, here::here("code", "report_gene.Rmd"), output_file)
}

#render_gene_report(input = "SST", type = "gene", output_file = "sst.zip")
#render_gene_report(input = "0060148", type = "pathway", output_file = "0060148.zip")
#render_gene_report(input = c("GSS", "SST"), type = "gene_list")


#logic to matching query type to rendered content
render_report_to_file <- function(data_values, file) {
  if (data_values$type == "gene" | data_values$type == "pathway" | data_values$type == "gene_list") {
    render_gene_report(data_values, output_file = file)
  } else {
    stop("no report for you!")
  }
}

#render_report_to_file(input = "GSS", type = "gene", file = "gss_trial.pdf")

