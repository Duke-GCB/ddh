library(here)
library(tidyverse)

source(here::here("code/generate_depmap_pathways.R"), local=TRUE)

test_that("enrichr_loop with empty gene_list can be arranged by Adjusted.P.value", {
  # test for error: "arrange() failed at implicit mutate() step.
  #   Could not create a temporary column for `Adjusted.P.value`."
  result = enrichr_loop(c(), focused_lib) %>%
    arrange(Adjusted.P.value)
  expect_equal(nrow(result), 0)
})
