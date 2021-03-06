#twitter_bot R script

library(tidyverse)
library(rtweet)
library(here)

#read current release information
source(here::here("code", "current_release.R"))

#load data
surprise_genes <- readRDS(file=here::here("data", paste0(release, "_surprise_genes.Rds")))
gene_summary <- readRDS(here::here("data", paste0(release, "_gene_summary.Rds")))
achilles <- readRDS(file=here::here("data", paste0(release, "_achilles.Rds")))
expression_join <- readRDS(file=here::here("data", paste0(release, "_expression_join.Rds")))
master_bottom_table <- readRDS(file=here::here("data", paste0(release, "_master_bottom_table.Rds")))
master_top_table <- readRDS(file=here::here("data", paste0(release, "_master_top_table.Rds")))

#load functions
source(here::here("code", "token.R"))
source(here::here("code", "fun_plots.R"))
source(here::here("code", "fun_graphs.R"))

twitter_save <- function(tmp_file, plot_id) {
  ggsave(tmp_file, plot_id, width = 16, height = 9, units = "cm", dpi = 150, device = "png")
}

#generate content
make_tweet <- function(gene_symbol) { #, random = TRUE
  gene_symbol <- sample(surprise_genes, 1)
  
  twitter_summary <- gene_summary %>% 
    filter(approved_symbol == gene_symbol)
  
  adverb_list <- c("an interesting", "an unexpected", "an unusual", "a surprising")
  adverb <- sample(adverb_list, 1)
  
  text <- paste0(twitter_summary$approved_symbol, ": ", twitter_summary$approved_name, " has ", adverb, " association with other genes. ")
  
  p1 <- make_celldeps(celldeps_data = achilles, expression_data = expression_join , gene_symbol, mean = 0)
  p1 <- p1 +
    labs(title = "Cell Line Dependency Curve", subtitle = "Each point shows a normalized genetic dependency score for a cell line")
  celldeps <- tempfile(fileext = ".png")
  twitter_save(celldeps, plot = p1)
  
  p2 <- make_cellbins(cellbins_data = achilles, expression_data = expression_join, gene_symbol)
  p2 <- p2 +
    labs(title = "Distribution of Dependency Scores", subtitle = "Kernel Density Estimate of the histogram of dependency scores")
  cellbins <- tempfile(fileext = ".png")
  twitter_save(cellbins, plot = p2)
  
  p3 <- make_lineage(celldeps_data = achilles, expression_data = expression_join, gene_symbol)
  p3 <- p3 +
    labs(title = "Cell Lineage Dependencies", subtitle = "Lineage dependency score summaries") +
    theme(axis.text.y = element_text(size = 8))
  lineage <- tempfile(fileext = ".png")
  twitter_save(lineage, plot = p3)
  
  p5 <- make_graph_report(toptable_data = master_top_table, bottomtable_data = master_bottom_table, gene_symbol, threshold = 10, deg = 2)
  p5 <- p5 +
    labs(title = "Network Graph of Dependencies", subtitle = "Shared top and bottom 10 dependency correlations")
  graph <- tempfile(fileext = ".png")
  twitter_save(graph, plot = p5)
  
  gene_symbol_url <- paste0("https://www.datadrivenhypothesis.org?show=detail&content=gene&symbol=", gene_symbol)
  
  text2 <- paste0("Check it out at ", gene_symbol_url)
  
  tweet_body <- str_c(text, text2) #cat(paste(text, text2, sep="\n"))
  tweet_list <- c(tweet_body, celldeps, cellbins, lineage, graph)
  
  return(tweet_list)
}

#post tweet
post_tweet(status = tweet_body, token = token, media = c(celldeps, cellbins, lineage, graph))


