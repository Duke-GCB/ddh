#load libraries
library(tidyverse)
library(here)
library(janitor)
library(corrr)
library(moderndive)
library(purrr)

#rm(list=ls()) 

#read current release information to set parameters for download
source(here::here("code", "current_release.R"))

time_begin_data <- Sys.time()

#These are log-fold change collapsed replicates with outliers and controls removed of PRISM
prism <- read_csv(prism_url, col_names = TRUE) %>% 
  clean_names()

#censor drugs that have gene names #"lta", "nppb", "prima1"
censor_names <- c("LTA", "NPPB", "PRIMA1")
censor_ids <- c("brd_k52914072_001_01_5_2_5_mts004", "brd_k89272762_001_12_7_2_5_hts", "brd_k15318909_001_10_5_2_5_hts")

prism <- 
  prism %>% 
  select(-any_of(censor_ids))

#get meta file
prism_meta <- read_csv(prismmeta_url, col_names = TRUE) %>% 
  clean_names() %>% 
  mutate(clean_drug = make_clean_names(column_name)) %>% 
  distinct(name, .keep_all = TRUE) %>%  #drop rows that have duplicate names
  filter(!name %in% censor_names) #censor 3 drugs from meta

#get CIDs into meta
url_exists <- function(x, non_2xx_return_value = FALSE, quiet = FALSE,...) {
  #from https://stackoverflow.com/questions/52911812/check-if-url-exists-in-r
  suppressPackageStartupMessages({
    require("httr", quietly = FALSE, warn.conflicts = FALSE)
  })
  
  sHEAD <- safely(httr::HEAD)
  sGET <- safely(httr::GET)
  
  # Try HEAD first since it's lightweight
  res <- sHEAD(x, ...)
  
  if (is.null(res$result) || 
      ((httr::status_code(res$result) %/% 200) != 1)) {
    
    res <- sGET(x, ...)
    
    if (is.null(res$result)) return(NA) # or whatever you want to return on "hard" errors
    
    if (((httr::status_code(res$result) %/% 200) != 1)) {
      if (!quiet) warning(sprintf("Requests for [%s] responded but without an HTTP status code in the 200-299 range", x))
      return(non_2xx_return_value)
    }
    
    return(TRUE)
    
  } else {
    return(TRUE)
  }
  
}
#from here:https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest$_Toc494865554
get_cid <- function(compound_name) {
  url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/", compound_name,"/property/IUPACName,MolecularFormula,CanonicalSMILES,ExactMass/CSV")
  if(url_exists(url) == TRUE) {
  tmp <- read_csv(url, col_types = cols()) #suppress msg
  tmp <- tmp %>% 
    mutate(name = compound_name) %>% 
    select(name, everything())
  } else {
    tmp <- 
      tibble(
        name = compound_name,
        CID = NA,
        IUPACName = NA,
        MolecularFormula = NA,
        CanonicalSMILES = NA,
        ExactMass = NA
      )
  }
return(tmp)
  }

cids <- tibble(
  name = character(),
  CID = double(),
  IUPACName = character(),
  MolecularFormula = character(),
  CanonicalSMILES = character(),
  ExactMass = double()
)

compounds <- prism_meta$name
#compounds <- sample(prism_meta$name, 10)

for (i in compounds) {
  tmp <- get_cid(i)
  cids <-
    cids %>%
    bind_rows(tmp)
  Sys.sleep(0.15) #per https://pubchemdocs.ncbi.nlm.nih.gov/programmatic-access
}

cids <- #some records return multiple; fortunately, best match is returned first
  cids %>% 
  distinct(name, .keep_all = TRUE) %>% 
  clean_names()

prism_meta <- 
  prism_meta %>% 
  left_join(cids, by = c("name" = "name")) %>%  #"smiles" = "CanonicalSMILES"
  filter(!is.na(name))

#make name/join/search df
prism_names <- prism_meta %>% 
  select("name", "moa", "cid", "clean_drug")

prism_long <- prism %>% #need this for joining below
  pivot_longer(cols = -x1, names_to = "drug", values_to = "log2fc") %>% 
  left_join(prism_names, by = c("drug" = "clean_drug")) %>% 
  filter(!is.na(name)) %>% #remove drugs which don't have a searchable name
  select(x1, name, log2fc)

#get log2FC values of ACHILLES for integration
achilles_log2fc_raw <- read_csv(achilles_log_url, col_names = TRUE)
achilles_guide_map <- read_csv(achilles_guide_map_url, col_names = TRUE)
achilles_rep_map <- read_csv(achilles_rep_map_url, col_names = TRUE)

#clean
achilles_guide_map$gene <- str_remove_all(achilles_guide_map$gene, "\\s\\(\\d+\\)")

achilles_log2fc <- achilles_guide_map %>% 
  left_join(achilles_log2fc_raw, by = c("sgrna" = "X1")) #%>% select(1:100)

achilles_log2fc_long <- achilles_log2fc %>% 
  pivot_longer(cols = "143B-311Cas9_RepA_p6_batch3":"BT549-311cas9 Rep A p5_batch2", names_to = "cell_line", values_to = "log2fc")

achilles_log2fc_long <- achilles_log2fc_long %>% 
  left_join(achilles_rep_map, by = c( "cell_line" = "replicate_ID")) %>% 
  select(DepMap_ID, gene, log2fc)

achilles_log2fc_long_mean <- achilles_log2fc_long %>% 
  group_by(DepMap_ID, gene) %>% 
  summarize(meanlog2fc = mean(log2fc), 
            sdlog2fc = sd(log2fc)) %>% 
  ungroup()

#combine, and pivot_wider for
combined <- achilles_log2fc_long_mean %>% 
  select(x1 = DepMap_ID, name = gene, log2fc = meanlog2fc) %>% 
  bind_rows(prism_long) %>% 
  rename(DepMap_ID = x1) %>% 
  pivot_wider(names_from = name, values_from = log2fc)

#Combined CORRELATION MATRIX
combined_cor <- combined %>%
  select(-DepMap_ID) %>% 
  correlate()

#test gene corrs ##comment out
# fav_gene <- "TSC1"
# combined_cor %>%
#   focus(!!fav_gene) %>%
#   arrange(desc(.[[2]]))

#resample for stats
#make some long files
combined_cor_long <- combined_cor %>% 
  stretch()

#Permutation tests
virtual_prism <- combined_cor_long %>% 
  filter(!is.na(r)) %>%   
  rep_sample_n(size = 20000, reps = 1000) %>%
  group_by(replicate) %>% 
  summarize(mean = mean(r), max = max(r), min = min(r), sd = sd(r))

mean_virtual_prism <- mean(virtual_prism$mean)
sd_virtual_prism <- mean(virtual_prism$sd)

drug_sd_threshold <- 2

prism_upper <- mean_virtual_prism + drug_sd_threshold*sd_virtual_prism
prism_lower <- mean_virtual_prism - drug_sd_threshold*sd_virtual_prism

#save
saveRDS(drug_sd_threshold, file = here::here("data", paste0(release, "_drug_sd_threshold.Rds")))
saveRDS(prism_lower, file = here::here("data", paste0(release, "_prism_lower.Rds")))
saveRDS(prism_upper, file = here::here("data", paste0(release, "_prism_upper.Rds")))
saveRDS(mean_virtual_prism, file = here::here("data", paste0(release, "_mean_virtual_prism.Rds")))
saveRDS(sd_virtual_prism, file = here::here("data", paste0(release, "_sd_virtual_prism.Rds")))

#cutoff and make tables
gene_drugs_table <- tibble(
  fav_gene = character(), 
  data = list()
)
drug_genes_table <- tibble(
  fav_drug = character(), 
  data = list()
)

#define list
#genes <- sample(names(select(combined_cor, A1BG:ZZZ3)), size = 10) #comment this out
#drugs <- sample(names(select(combined_cor, !rowname:ZZZ3)), size = 10) #comment this out
genes <- names(select(combined_cor, A1BG:ZZZ3))
drugs <- names(select(combined_cor, !rowname:ZZZ3))
#18524+4514 == 23038 (same as combined_cor)

#drug table for a gene
for (fav_gene in genes) {
  message("Drug tables for ", fav_gene)
  gene_top <- 
    combined_cor %>% 
    focus(fav_gene) %>% 
    arrange(desc(.[[2]])) %>% #use column index
    filter(rowname %in% drugs, #remove genes
           .[[2]] > prism_upper) %>% #mean +/- 2sd
    rename(drug = rowname, 
         r2 = fav_gene) %>% 
    mutate(r2 = round(r2, 2), 
           z_score = round((r2 - mean_virtual_prism)/sd_virtual_prism, 1)) %>% 
    select(drug, z_score, r2)

  gene_table <- 
    gene_top %>% 
    mutate(fav_gene = fav_gene) %>% 
    group_by(fav_gene) %>% 
    nest()
  
  gene_drugs_table <- gene_drugs_table %>% 
    bind_rows(gene_table)
}

#gene table for a drug query
for (fav_drug in drugs) {
  message("Gene tables for ", fav_drug)
  drug_top <- 
    combined_cor %>% 
    focus(fav_drug) %>% 
    arrange(desc(.[[2]])) %>% #use column index
    filter(rowname %in% genes, #remove drugs
           .[[2]] > prism_upper) %>% #mean +/- 2sd
    rename(gene = rowname, 
           r2 = fav_drug) %>% 
    mutate(r2 = round(r2, 2), 
           z_score = round((r2 - mean_virtual_prism)/sd_virtual_prism, 1)) %>% 
    select(gene, z_score, r2)
  
  drug_table <- drug_top %>% 
    mutate(fav_drug = fav_drug) %>% 
    group_by(fav_drug) %>% 
    nest()
  
  drug_genes_table <- drug_genes_table %>% 
    bind_rows(drug_table)
}
  
#TEST get data out
# make_drug_table <- function(gene_data = gene_drugs_table, gene_symbol) {
#   gene_data %>%
#     dplyr::filter(fav_gene %in% gene_symbol) %>%
#     tidyr::unnest(data) %>%
#     dplyr::arrange(desc(r2)) %>%
#     dplyr::rename("Query" = "fav_gene", "Drug" = "drug", "R^2" = "r2", "Z Score" = "z_score")
# }
# make_gene_table <- function(drug_data = drug_genes_table, drug_name) {
#   drug_data %>%
#     dplyr::filter(fav_drug %in% drug_name) %>%
#     tidyr::unnest(data) %>%
#     dplyr::arrange(desc(r2)) %>%
#     dplyr::rename("Query" = "fav_drug", "Gene" = "gene", "R^2" = "r2", "Z Score" = "z_score")
# }
#combined_cor_long %>% arrange(desc(r)) %>% filter(x %in% drugs) %>% filter(y %in% genes)

#save files
saveRDS(prism, file = here::here("data", paste0(release, "_prism.Rds")))
saveRDS(prism_meta, file = here::here("data", paste0(release, "_prism_meta.Rds")))
saveRDS(prism_names, file = here::here("data", paste0(release, "_prism_names.Rds")))
#saveRDS(combined_cor, file = here::here("data", paste0(release, "_combined_cor.Rds")))
saveRDS(gene_drugs_table, file=here::here("data", paste0(release, "_gene_drugs_table.Rds")))
saveRDS(drug_genes_table, file=here::here("data", paste0(release, "_drug_genes_table.Rds")))


#how long
time_end_data <- Sys.time()
