library(tidyverse)
library(cowplot)
library(scico)
library(plotly)
library(gganatogram)
library(ragg)

make_cellbins <- function(cellbins_data = achilles, expression_data = expression_join, gene_symbol) {
  plot_complete <- 
    cellbins_data %>% #plot setup
    select(X1, any_of(gene_symbol)) %>%
    left_join(expression_data, by = "X1") %>%
    select(-X1) %>%
    pivot_longer(cols = where(is.numeric), names_to = "gene_symbol", values_to = "dep_score") %>% 
    group_by(gene_symbol) %>% 
    arrange(dep_score) %>% 
    mutate(med = median(dep_score, na.rm = T)) %>% 
    ungroup() %>% 
    filter(!is.na(dep_score)) %>% 
    ggplot() +
    geom_vline(xintercept = 1, color = "lightgray") +
    geom_vline(xintercept = -1, color = "lightgray") +
    geom_vline(xintercept = 0) +
    geom_linerange(aes(xmin = -Inf, xmax = med, 
                       y = fct_reorder(gene_symbol, med), 
                       color = fct_reorder(gene_symbol, med)),
                   linetype = "dotted",
                   size = .2) +
    ggdist::stat_halfeye(aes(x = dep_score, y = fct_reorder(gene_symbol, med),
                             color = fct_reorder(gene_symbol, med), 
                             fill = after_scale(colorspace::lighten(color, .6, space = "HLS")),
                             point_fill = after_scale(colorspace::lighten(color, .5, space = "HLS"))),
                         .width = c(.5, .95),
                         shape = 21,
                         stroke = .7,
                         point_size = 2) +
    geom_vline(xintercept = 0, alpha = .2) +
    labs(x = "Dependency Score (binned)", y = NULL, color = "Query \nGene", fill = "Query \nGene") +
    scale_y_discrete(expand = c(.03, .03)) +
    scale_color_scico_d(palette = "lapaz", guide = "legend", end = .8) +
    scale_fill_scico_d(palette = "lapaz", guide = "legend", end = .8) +
    scale_fill_viridis(discrete = TRUE, option = "D", direction = 1, guide = "legend") +
    guides(
      color = guide_legend(size = 1, reverse = T),
      fill = guide_legend(size = 1, reverse = T)
    ) +
    theme_cowplot(font_size = 16) +
    theme(legend.position = "none", axis.line.y = element_blank(), axis.ticks.y = element_blank(), 
          axis.text.y = element_text(size = 17), text = element_text(family = "Chivo")) +
    NULL
  
  if(length(gene_symbol) == 1){
    plot_complete  <- plot_complete +
      guides(fill = "none")
  } else {
    plot_complete
  }
  return(plot_complete)
}

#figure legend
plot_cellbins_title <- "Kernel density estimate."
plot_cellbins_legend <- "A smoothed version of the histogram of Dependency Scores. Dependency scores across all cell lines for queried genes, revealing overall influence of a gene on cellular fitness"

make_celldeps <- function(celldeps_data = achilles, expression_data = expression_join, gene_symbol, mean) {
  plot_complete <- 
    celldeps_data %>% #plot setup
    select(X1, any_of(gene_symbol)) %>%
    left_join(expression_data, by = "X1") %>%
    select(-X1) %>%
    pivot_longer(cols = where(is.numeric), names_to = "gene_symbol", values_to = "dep_score") %>% 
    group_by(gene_symbol) %>% 
    arrange(dep_score) %>% 
    mutate(
      rank = 1:n(),
      med = median(dep_score, na.rm = T)
    ) %>% 
    ungroup() %>% 
    ggplot(aes(x = rank, 
               y = dep_score, 
               text = paste0("Cell Line: ", cell_line), 
               color = fct_reorder(gene_symbol, med),
               fill = fct_reorder(gene_symbol, med)
      )) +
      geom_hline(yintercept = mean) +
      geom_hline(yintercept = 1, color = "lightgray", linetype = "dashed") +
      geom_hline(yintercept = -1, color = "lightgray", linetype = "dashed") +
      geom_hline(yintercept = 0) +
      geom_point(size = 1, stroke = .1, alpha = 0.4) + 
      scale_x_discrete(expand = expansion(mult = 0.02), na.translate = FALSE) +
      scale_color_scico_d(palette = "lapaz", guide = "legend", end = .8) +
      scale_fill_scico_d(palette = "lapaz", guide = "legend", end = .8) +
      guides(
        color = guide_legend(reverse = T, override.aes = list(size = 5)),
        fill = guide_legend(reverse = T, override.aes = list(size = 5))
      ) +
      labs(x = "Rank", y = "Dependency Score", color = "Query \nGene", fill = "Query \nGene") +
      theme_cowplot(font_size = 16) +
      theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.line.x = element_blank(), text = element_text(family = "Chivo")) + # axis.title.x=element_blank()
      NULL
  
  
  if(length(gene_symbol) == 1){
    plot_complete  <- plot_complete +
      guides(color = "none")
  } else {
    plot_complete
  }
  return(plot_complete)
}

#figure legend
plot_celldeps_title <- "Cell Line Dependency Curve."
plot_celldeps_legend <- "Each point shows the ranked dependency score for a given cell line. Cells with dependency scores less than -1 indicate a cell that the query gene is essential within. Cells with dependency scores close to 0 show no changes in fitness when the query gene is knocked out. Cells with dependency scores greater than 1 have a gain in fitness when the query gene is knocked-out"

# make cell anatogram
make_cellanatogram <- function(cellanatogram_data = subcell, gene_symbol) {
  plot_complete <- cellanatogram_data %>% 
    filter_all(any_vars(gene_name %in% gene_symbol)) %>% 
    filter(!is.na(type)) %>% 
    select(-value) %>% 
    add_count(main_location) %>% 
    mutate(value = as_factor(n)) %>% 
    gganatogram(outline = T, fillOutline='grey95', organism="cell", fill = "value") +
    theme_void() +  
    coord_fixed() +
    scale_fill_viridis(discrete = TRUE) +
    labs(fill = "Count") +
    theme(text = element_text(family = "Chivo"))
  
  if(length(gene_symbol) == 1){
    plot_complete  <- plot_complete +
      guides(fill = "none")
  } else {
    plot_complete
  }
  return(plot_complete)
}

# make lineage plot
make_lineage <- function(celldeps_data = achilles, expression_data = expression_join, gene_symbol) {
  data_full <- celldeps_data %>% #plot setup
    select(X1, any_of(gene_symbol)) %>%
    left_join(expression_data, by = "X1") %>%
    select(-X1) %>%
    pivot_longer(cols = where(is.numeric), names_to = "gene_symbol", values_to = "dep_score") %>%
    dplyr::mutate_at("lineage", function(str) {
      str <- str_replace_all(str, "\\_", " ")
      str <- str_to_title(str)
      return(str)
    }) %>% 
    drop_na(lineage) %>% 
    drop_na(dep_score) %>% 
    group_by(lineage) %>% 
    mutate(mean = mean(dep_score)) %>% 
    ungroup() %>% 
    mutate(lineage = fct_reorder(lineage, mean))
  
  data_mean <- data_full %>% 
    group_by(lineage) %>% 
    summarize(dep_score = mean(dep_score))
  
  plot_complete <- 
    ggplot() +
      geom_vline(xintercept = 0) +
      geom_linerange(data = data_mean,
                     aes(xmin = -Inf, xmax = dep_score, y = lineage),
                     color = "grey60",
                     linetype = "dotted") +
      ggdist::stat_interval(data = data_full,
                            aes(x = dep_score, y = lineage),
                            orientation = "horizontal",
                           .width = c(.05, .5, .95)
      ) +
      geom_point(data = data_mean, aes(x = dep_score, y = lineage), color = "grey20") +
      scale_color_manual(values = c("#aae3dd", "#19acb5", "#036d77"), labels = c("95%", "50%", "5%"), name = "") +
      guides(color = guide_legend(reverse = TRUE)) +
      labs(x = "Dependency Score", y = NULL, title = "Cell Lineage:") +
      theme_cowplot(font_size = 16) +
      theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), 
            text = element_text(family = "Chivo"), legend.position = "bottom", 
            plot.title = element_text(size = 14), plot.title.position = "plot")
  return(plot_complete)
}

#figure legend
plot_celllin_title <- "Cell Line Lineage Dependencies"
plot_celllin_legend <- "Each point shows the mean dependency score for a given cell lineage, with box plots showing median, interquartile ranges, and outliers."

# make sublineage plot
make_sublineage <- function(celldeps_data = achilles, expression_data = expression_join, gene_symbol) {
  data_full <- celldeps_data %>% #plot setup
    select(X1, any_of(gene_symbol)) %>%
    left_join(expression_data, by = "X1") %>%
    select(-X1) %>%
    pivot_longer(cols = where(is.numeric), names_to = "gene_symbol", values_to = "dep_score") %>% 
    dplyr::mutate_at("lineage_subtype", function(str) {
      str <- str_replace_all(str, "\\_", " ")
      str <- if_else(str_detect(str, "^[:lower:]"), str_to_title(str), str)
      return(str)
    })  %>% 
    drop_na(lineage_subtype) %>% 
    drop_na(dep_score) %>% 
    group_by(lineage_subtype) %>% 
    mutate(mean = mean(dep_score)) %>% 
    ungroup() %>% 
    mutate(lineage_subtype = fct_reorder(lineage_subtype, mean))
  
  data_mean <- data_full %>% 
    group_by(lineage_subtype) %>% 
    summarize(dep_score = mean(dep_score))
  
  plot_complete <- 
    ggplot() +
    geom_vline(xintercept = 0) +
    geom_linerange(data = data_mean,
                   aes(xmin = -Inf, xmax = dep_score, y = lineage_subtype),
                   color = "grey60",
                   linetype = "dotted") +
    ggdist::stat_interval(data = data_full,
                          aes(x = dep_score, y = lineage_subtype),
                          orientation = "horizontal",
                          .width = c(.05, .5, .95)
    ) +
    geom_point(data = data_mean, aes(x = dep_score, y = lineage_subtype), color = "grey20") +
    scale_color_manual(values = c("#aae3dd", "#19acb5", "#036d77"), labels = c("95%", "50%", "5%"), name = "") +
    guides(color = guide_legend(reverse = TRUE)) +
    labs(x = "Dependency Score", y = NULL, title = "Cell Sublineage:") +
    theme_cowplot(font_size = 16) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), 
          text = element_text(family = "Chivo"), legend.position = "bottom", 
          plot.title = element_text(size = 14), plot.title.position = "plot")
  return(plot_complete)
}

#figure legend
plot_cellsublin_title <- "Cell Line Sub-Lineage Dependencies"
plot_cellsublin_legend <- "Each point shows the mean dependency score for a given cell sublineage, with box plots showing median, interquartile ranges, and outliers."

