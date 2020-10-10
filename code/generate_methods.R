#upon each release, run

#read current release information 
source(here::here("code", "current_release.R"))

source(here::here("code", "find_threshold.R"))
source(here::here("code", "generate_depmap_data.R"))
source(here::here("code", "generate_pubmed_data.R"))
source(here::here("code", "generate_depmap_stats.R"))

#generate table data
#go to generate_depmap_tables.R
#set methods = TRUE in header, and then source

#generate pathway data
source(here::here("code", "generate_depmap_pathways.R"))
master_positive <- generate_positive_data(gene_group = c("TP53"), achilles_cor = achilles_cor, achilles_upper = achilles_upper, gene_summary = gene_summary)
saveRDS(master_positive, file=here::here("data", paste0(release, "_", master_positive_filename)))
master_negative <- generate_negative_data(gene_group = c("TP53"), achilles_cor = achilles_cor, achilles_lower = achilles_lower, gene_summary = gene_summary)
saveRDS(master_negative, file=here::here("data", paste0(release, "_", master_negative_filename)))

#knit methods