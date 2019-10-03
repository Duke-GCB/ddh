library(tidyverse)
library(lubridate)
library(here)
library(janitor)
library(feather)
library(corrr)
library(purrr)

params <- list(
  release="19Q3",
  datadir="data",
  achilles_url="https://ndownloader.figshare.com/files/16757666",
  ccle_url="https://ndownloader.figshare.com/files/16757690",
  cclemeta_url="https://ndownloader.figshare.com/files/16757723",
  na_cutoff=400
)

download_achilles <- function(params) {
  read_csv(params$achilles_url, col_names = TRUE) %>%
      `colnames<-`(str_remove_all(names(.), "\\s\\(\\d+\\)"))
}

download_expression <- function(params) {
  read_csv(params$ccle_url, col_names = TRUE) %>%
    `colnames<-`(str_remove_all(names(.), "\\s\\(\\d+\\)"))
}

download_expression_metadata <- function(params) {
  read_csv(params$cclemeta_url, col_names = TRUE) %>%
    clean_names
}

create_achilles_long <- function(achilles) {
  achilles %>%
    gather("gene", "dep_score", -X1)
}

create_expression_long <- function(expression, achilles) {
  expression %>%
    filter(expression$X1 %in% achilles$X1 == TRUE) %>% #matches cells
    gather("gene", "gene_expression", -X1) %>%
    arrange(desc(gene_expression))
}

create_achilles_cor <- function(params, expression, achilles) {
  achilles_long <- create_achilles_long(achilles)
  expression_long <- create_expression_long(expression, achilles)

  no_expression <- expression_long %>%
    filter(gene_expression == 0) %>%
    unite(X1, gene, col = "match", sep = "-", remove = TRUE) %>%
    pull(match)

  #make new match df
  achilles_no0 <- achilles_long %>%
    unite(X1, gene, col = "match", sep = "-", remove = FALSE) %>%
    filter(match %in% no_expression == FALSE) %>%
    select(-match) %>%
    spread(gene, dep_score)

  na_cutoff <- params$na_cutoff

  toomanyNAs <- achilles_no0 %>%
    summarise_all(list(~sum(is.na(.)))) %>%
    gather(gene, NAs) %>%
    arrange(desc(NAs)) %>%
    filter(NAs > na_cutoff) %>%
    pull(gene)

  achilles_clean <- achilles_no0 %>%
    select(-one_of(toomanyNAs)) #4491 elements for NAs > 100; 1871 elements for NAs > 400

  achilles_clean %>% #originally 'achilles'
    select(-X1) %>%
    correlate() #(diagonal = 0) set to 0 so easy to summarize, but should be NA; so added na.rm = TRUE to fun() in EDA
}

create_achilles_cor_long <- function(achilles_cor) {
  achilles_cor %>%
    stretch() #310M observations across 3 variables (x, y, r)
}

create_expression_cor <- function(expression) {
  expression %>%
    select(-X1) %>%
    correlate(diagonal = 0) #set to 0 so easy to summarize
}

create_expression_cor_long <- function(expression_cor) {
  expression_cor %>%
    stretch() #310M observations across 3 variables (x, y, r)
}

save_dataframe <- function(params, df, filename_suffix) {
  write_feather(df, path = here(params$datadir, paste0(params$release, filename_suffix)))
}

correlate_and_save_achilles <- function (params, achilles, expression) {
  achilles_cor <- create_achilles_cor(params, expression, achilles)
  save_dataframe(params, achilles_cor,  "_achilles_cor.feather")
}

correlate_and_save_expression <- function (params, expression) {
  expression_cor <- create_expression_cor(expression)
  save_dataframe(params, expression_cor, "_expression_cor.feather")
}

correlate_and_save <- function(params) {
  achilles <- download_achilles(params)
  save_dataframe(params, achilles, "_achilles.feather")

  expression <- download_expression(params)
  save_dataframe(params, expression, "_expression.feather")

  expression_id <- download_expression_metadata(params)
  save_dataframe(params, expression_id, "_expression_id.feather")

  correlate_and_save_achilles(params, expression, achilles)
}

start_time <- Sys.time()

print(correlate_and_save(params))

end_time <- Sys.time()
time_taken <- round(as.duration(start_time %--% end_time)/dminutes(1), digits = 1)
print(time_taken)
