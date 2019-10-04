library(tidyverse)
library(lubridate)
library(here)
library(readxl)
library(janitor)
library(feather)
library(corrr)
library(purrr)
library(viridis)
library(enrichR)
library(glue)
library(moderndive)
library(rmarkdown)
library(gridExtra)
library(knitr)
library(pander)
library(vroom)

params <- list(
  release="19Q3",
  datadir="data",
  gene_name="BRCA1"
)

focused_lib = readLines(here::here(params$datadir, 'focused_lib.csv'))

enrichr_loop <- function(gene_list, databases){
  if(is_empty(gene_list)){
    flat_complete <- NULL
    return(flat_complete)
  } else {
    flat_complete <- as_tibble()
    for(lib in databases){
      enriched <- enrichr(gene_list, lib)

      flat <- flatten_dfc(enriched) %>%
        mutate(enrichr = lib)

      flat_complete <- flat_complete %>%
        bind_rows(flat)
    }
    flat_complete <- flat_complete %>%
      arrange(Adjusted.P.value) %>%
      select(enrichr, Term, Overlap)

    flat_complete$enrichr <- str_replace_all(flat_complete$enrichr, "\\_", " ")
    flat_complete$Term <- str_replace_all(flat_complete$Term, "\\_", " ")
    return(flat_complete)
  }
}

load_dataframe <- function(params, filename_suffix) {
  read_feather(here::here(params$datadir, paste0(params$release, filename_suffix)))
}

run <- function(params) {
  achilles <- load_dataframe(params, "_achilles.feather")
  achilles_cor <- load_dataframe(params, "_achilles_cor.feather")
  class(achilles_cor) <- c("cor_df", "tbl_df", "tbl", "data.frame") #define class so functions (eg focus) can work on reloaded df
  achilles_cor_long <- stretch(achilles_cor)

  expression_id <- load_dataframe(params, "_expression_id.feather")
  expression_join <- expression_id %>%
    rename(X1 = dep_map_id) %>%
    select(X1, stripped_cell_line_name, lineage)

  virtual_achilles <- achilles_cor_long %>% #achilles_cor_long already has all of the variables in a long format
    filter(!is.na(r)) %>%
    rep_sample_n(size = 20000, reps = 1000) %>% #larger sample size, less error (but only 625 sets, and we're mimicking 1000?, but 310M combinations, so probably OK)
    group_by(replicate) %>%
    summarize(mean = mean(r), max = max(r), min = min(r), sd = sd(r)) #how to handle + vs. - correlation?

  mean_virtual_achilles <- mean(virtual_achilles$mean)
  sd_virtual_achilles <- mean(virtual_achilles$sd)

  sd_threshold <- 3

  achilles_upper <- mean_virtual_achilles + sd_threshold*sd_virtual_achilles
  achilles_lower <- mean_virtual_achilles - sd_threshold*sd_virtual_achilles

  id <- read_feather(here::here("data", "id.feather"))

  proteins <- id %>%
    select(gene, protein_name)

  #gene_group <- c(params$gene_name)
  achilles_names <- names(achilles)
  gene_group = achilles_names[achilles_names != "X1"]

  achilles_long <- achilles %>%
    gather("gene", "dep_score", -X1)

  for (fav_gene in gene_group) {
    if(fav_gene %in% names(achilles_cor) == 1){ #this code checks to see if the gene is in the analysis, and if not, skips
      dep_top <- achilles_cor %>%
        focus(fav_gene) %>%
        arrange(desc(.[[2]])) %>% #use column index
        filter(.[[2]] > achilles_upper) %>% #formerly top_n(20), but changed to mean +/- 3sd
        rename(gene = rowname) %>%
        left_join(proteins, by = "gene") %>%
        select(gene, protein_name, fav_gene) %>%
        rename(protein = protein_name, r2 = fav_gene)

      dep_bottom <- achilles_cor %>%
        focus(fav_gene) %>%
        arrange(.[[2]]) %>% #use column index
        filter(.[[2]] < achilles_lower) %>% #formerly top_n(20), but changed to mean +/- 3sd
        rename(gene = rowname) %>%
        left_join(proteins, by = "gene") %>%
        select("gene", "protein_name", fav_gene) %>%
        rename(protein = protein_name, r2 = fav_gene)
      #this is to get neg correlators

      #pathway enrichment analyses
      flat_top_complete <- dep_top %>%
        pull("gene") %>%
        c(fav_gene, .) %>%
        enrichr_loop(., focused_lib)

      #bottom
      flat_bottom_complete <- dep_bottom %>%
        pull("gene") %>%
        enrichr_loop(., focused_lib)

      #plot setup
      target_achilles <- achilles_long %>%
        filter(gene == fav_gene) %>%
        left_join(expression_join, by = "X1") %>%
        select(stripped_cell_line_name, lineage, dep_score)

      target_achilles_top <- target_achilles %>%
        top_frac(dep_score, n = 0.01)

      target_achilles_bottom <- target_achilles %>%
        top_frac(dep_score, n = -0.01) %>%
        arrange(dep_score)

      #plot1
      dep_plot1 <- ggplot(target_achilles) +
        geom_histogram(aes(x = dep_score), binwidth = 0.25, color = "lightgray") +
        labs(x = "Dependency Score (binned)") +
        theme_light()

      #plot2
      dep_plot2 <- ggplot(target_achilles) +
        geom_point(aes(x = fct_rev(fct_reorder(target_achilles$stripped_cell_line_name, target_achilles$dep_score, .desc = TRUE)), y = dep_score)) +
        labs(x = "Cell Lines", y = "Dependency Score") +
        geom_hline(yintercept = mean_virtual_achilles) +
        geom_hline(yintercept = achilles_upper, linetype="dashed") +
        geom_hline(yintercept = achilles_lower, linetype="dashed") +
        geom_hline(yintercept = 0) +
        geom_point(data = target_achilles_top, aes(x = stripped_cell_line_name, y = dep_score), color = "red") +
        geom_point(data = target_achilles_bottom, aes(x = stripped_cell_line_name, y = dep_score), color = "red") +
        theme_light() +
        theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + # axis.title.x=element_blank()
        NULL

      #summary
      fav_gene_entrez <- pull(id[match(fav_gene, id$gene), 7])
      if(is.na(fav_gene_entrez) == TRUE){
        lookup <- c("NULL")
        summary <- as_tibble(colnames("X1")) #make tibble to avoid atomic vector error in report
      } else {
        # paste into url
        lookup <- paste0("http://mygene.info/v3/gene/", fav_gene_entrez, "?fields=summary") #My gene info lookup; https://docs.mygene.info/en/latest/doc/annotation_service.html
        summary <- vroom(lookup, col_names = FALSE) %>% filter(str_detect(X1, "summary"))
      }

      #render output
      render("code/report_depmap_complete.Rmd", output_dir = here::here("results"), output_file = paste0(fav_gene, '_depmap.pdf'))
    } else {
      #summary
      fav_gene_entrez <- pull(id[match(fav_gene, id$gene), 7])
      if(is.na(fav_gene_entrez) == TRUE){
        lookup <- c("NULL")
        summary <- as_tibble(colnames("X1")) #make tibble to avoid atomic vector error in report
      } else {
        # paste into url
        lookup <- paste0("http://mygene.info/v3/gene/", fav_gene_entrez, "?fields=summary") #My gene info lookup; https://docs.mygene.info/en/latest/doc/annotation_service.html
        summary <- vroom(lookup, col_names = FALSE) %>% filter(str_detect(X1, "summary"))
      }

      #render output
      render("code/report_dummy_depmap.Rmd", output_dir = here::here("results"), output_file = paste0(fav_gene, '_depmap.pdf'))
    }
  }
}


run(params)
