library(tidyverse)
library(tidygraph)
library(networkD3)
library(ggraph)
library(viridis)
library(visNetwork)

source(here::here("code", "fun_tables.R")) #for make_table funs

setup_graph <- function(toptable_data = master_top_table, bottomtable_data = master_bottom_table, gene_symbol, threshold = 10) {
  #make empty tibble
  dep_network <- tibble()
  #either find top/bottom correlated genes if given single gene, or take list to fill gene_list
  if(length(gene_symbol) == 1){
    #find top and bottom correlations for fav_gene
    dep_top <- make_top_table(toptable_data, gene_symbol) %>%
      slice(1:threshold)
    
    dep_bottom <- make_bottom_table(bottomtable_data, gene_symbol) %>%
      slice(1:threshold) #limit for visualization?
    
    #this takes the genes from the top and bottom, and pulls them to feed them into a for loop
    gene_list <- dep_top %>%
      bind_rows(dep_bottom) %>%
      dplyr::pull("Gene")
  } else {
    gene_list <- toptable_data %>% #this code ensures that the list of genes from a pathway are in the data
      dplyr::filter(fav_gene %in% gene_symbol) %>%
      dplyr::pull(fav_gene)
  }
  #this loop will take each gene, and get their top and bottom correlations, and build a df containing the top n number of genes for each gene
  for (i in gene_list){
    message("Getting correlations from ", i)
    dep_top_related <- toptable_data %>%
      dplyr::filter(fav_gene == i) %>%
      tidyr::unnest(data) %>%
      dplyr::ungroup(.) %>% 
      dplyr::arrange(desc(r2)) %>%
      dplyr::slice(1:threshold) %>% 
      dplyr::mutate(x = i, origin = "pos") %>% 
      dplyr::rename(y = gene) %>%
      dplyr::select(x, y, r2, origin)
    
    dep_bottom_related <- bottomtable_data %>%
      dplyr::filter(fav_gene == i) %>%
      tidyr::unnest(data) %>%
      dplyr::ungroup(.) %>% 
      dplyr::arrange(r2) %>%
      dplyr::slice(1:threshold) %>% 
      dplyr::mutate(x = i, origin = "neg") %>% 
      dplyr::rename(y = gene) %>%
      dplyr::select(x, y, r2, origin)
    
    #each temp object is bound together, and then bound to the final df for graphing
    dep_related <- dep_top_related %>%
      bind_rows(dep_bottom_related)
    
    dep_network <- dep_network %>%
      bind_rows(dep_related)
  }
  return(dep_network)
}
#tests
#setup_graph(gene_symbol = "SDHA")
#setup_graph(gene_symbol = c("SDHA", "SDHB"))
#setup_graph(gene_symbol = c("GSS", "SST"))

make_graph <- function(toptable_data = master_top_table, bottomtable_data = master_bottom_table, gene_symbol, threshold = 10, deg = 2) {
  #get dep_network object
  dep_network <- setup_graph(toptable_data, bottomtable_data, gene_symbol, threshold)
  
  if(length(gene_symbol) == 1){
    dep_top <- make_top_table(toptable_data, gene_symbol) %>% slice(1:threshold) #redundant with above, but need these objs
    dep_bottom <- make_bottom_table(bottomtable_data, gene_symbol) %>% slice(1:threshold)
  } else {
    dep_network_top <- dep_network %>% filter(origin == "pos") %>% pull(y)
    dep_network_bottom <- dep_network %>% filter(origin == "neg") %>% pull(y)
  }
  
  #make graph
  graph_network <- tidygraph::as_tbl_graph(dep_network)
  if(length(gene_symbol) == 1){
    nodes <-  as_tibble(graph_network) %>%
      rowid_to_column("id") %>%
      mutate(degree = igraph::degree(graph_network),
             group = case_when(name %in% gene_symbol == TRUE ~ "Query Gene", 
                               name %in% dep_top$Gene == TRUE ~ "Positive",
                               name %in% dep_bottom$Gene == TRUE ~ "Negative",
                               TRUE ~ "Connected"), 
             group = as_factor(group), 
             group = fct_relevel(group, c("Query Gene", "Positive", "Negative", "Connected")))  %>%
      arrange(group)
  } else {
    nodes <-  as_tibble(graph_network) %>%
      rowid_to_column("id") %>%
      mutate(degree = igraph::degree(graph_network),
             group = dplyr::case_when(name %in% gene_symbol == TRUE ~ "Query Gene", 
                                      name %in% dep_network_top == TRUE ~ "Positive",
                                      name %in% dep_network_bottom == TRUE ~ "Negative"),
             group = as_factor(group), 
             group = fct_relevel(group, c("Query Gene", "Positive", "Negative")))  %>% #you don't end up with "connected" in a multi-gene list
      arrange(group) 
  }
  
  links <- graph_network %>%
    activate(edges) %>% # %E>%
    as_tibble()
  
  # determine the nodes that have at least the minimum degree
  nodes_filtered <- nodes %>%
    filter(degree >= deg) %>%  #input$degree
    as.data.frame
  
  # filter the edge list to contain only links to or from the nodes that have the minimum or more degree
  links_filtered <- links %>%
    filter(to %in% nodes_filtered$id & from %in% nodes_filtered$id) %>%
    as.data.frame
  
  # re-adjust the from and to values to reflect the new positions of nodes in the filtered nodes list
  links_filtered$from <- match(links_filtered$from, nodes_filtered$id) - 1
  links_filtered$to <- match(links_filtered$to, nodes_filtered$id) - 1
  
  #check to see if setting degree removed all links; if so, then throws error, so this fills a dummy links_filtered df to plot only nodes
  if(nrow(links_filtered) == 0) {links_filtered <- tibble("from" = 0, "to" = 0, "r2" = 1, "origin" = "pos")}
  
  #use color meter to get hexdec color values
  node_color <- 'd3.scaleOrdinal(["#EDA555","#AD677D", "#0C2332", "#544097"])'
  
  #check to see if query gene is missing; if so, then adds a dummy so it shows up on graph, but disconnected
  if(sum(str_detect(nodes_filtered$group, "Query Gene")) == 0){
    dummy <- tibble("id" = max(nodes_filtered$id) + 1, "name" = gene_symbol, "degree" = 1, "group" = "Query Gene")
    nodes_filtered <- bind_rows(nodes_filtered, dummy)
    node_color <- 'd3.scaleOrdinal(["#AD677D", "#0C2332", "#544097", "#EDA555"])' #change color order for consistency
  }
  
  # Link: https://datadrivenhypothesis.com/?show=gene&query_type=gene&symbol=TSC1
  # Javascript for the tooltip upon clicking of a node
  # MyClickScript <- 'alert("You clicked " + d.name + " which is in row " +
  #      (d.index + 1) +  " of your original R data frame");'
  MyClickScript <- 'window.open("https://datadrivenhypothesis.com/?show=gene&query_type=gene&symbol=" + d.name)'
  
  
  forceNetwork(Links = links_filtered, 
               Nodes = nodes_filtered, 
               Source = "from", 
               Target ="to", 
               NodeID = "name", 
               Group = "group", 
               zoom = TRUE, 
               bounded = TRUE, 
               opacity = 0.8,
               opacityNoHover = 100, 
               Nodesize = "degree", 
               colourScale = node_color, 
               legend = TRUE,
               clickAction=MyClickScript)
}

# make_graph(gene_symbol = "SDHA")
# make_graph(gene_symbol = "AADACL4")
#make_graph(gene_symbol = c("SDHA", "SDHB"))
#make_graph(gene_symbol = c("GSS", "SST"))

#' Create network graph visualization using visNetwork
#' 
#' This function takes in dependency correlations and a gene query list to then output a dependency network graph 
#' visualization containing the top/bottom threshold for each of the top/bottom threshold of the gene query list 
#' using visNetwork.
#'
#' @param toptable_data A tibble of genes and their associated top correlated genes
#' @param bottomtable_data A tibble of genes and their associated bottom correlated genes
#' @param gene_symbol A character or character vector of gene_symbols used to create network graph
#' @param threshold A numerical representing the number of genes to pull from top and bottom tables
#' @param deg A numerical representing the minimum number of connections for a gene to be connected to the network
#'
#' @return NULL - Outputs a complete network graph
#' @export
#'
#' @examples
make_graph_visNetwork <- function(toptable_data = master_top_table, bottomtable_data = master_bottom_table, gene_symbol, threshold = 10, deg = 2) {
  #get dep_network object
  dep_network <- setup_graph(toptable_data, bottomtable_data, gene_symbol, threshold)
  
  # TODO: Output dep top objects from the setup_graph as well. Need to create an object in setup graph with each of these objs in it to access it.
  if(length(gene_symbol) == 1){
    dep_top <- make_top_table(toptable_data, gene_symbol) %>% slice(1:threshold) #redundant with above, but need these objs
    dep_bottom <- make_bottom_table(bottomtable_data, gene_symbol) %>% slice(1:threshold)
  } else {
    dep_network_top <- dep_network %>% filter(origin == "pos") %>% pull(y)
    dep_network_bottom <- dep_network %>% filter(origin == "neg") %>% pull(y)
  }
  
  #make graph
  graph_network <- tidygraph::as_tbl_graph(dep_network)
  if(length(gene_symbol) == 1){
    nodes <-  as_tibble(graph_network) %>%
      rowid_to_column("id") %>%
      mutate(degree = igraph::degree(graph_network),
             group = case_when(name %in% gene_symbol == TRUE ~ "Query Gene", 
                               name %in% dep_top$Gene == TRUE ~ "Positive",
                               name %in% dep_bottom$Gene == TRUE ~ "Negative",
                               TRUE ~ "Connected"), 
             group = as_factor(group), 
             group = fct_relevel(group, c("Query Gene", "Positive", "Negative", "Connected")))  %>%
      arrange(group)
  } else {
    nodes <-  as_tibble(graph_network) %>%
      rowid_to_column("id") %>%
      mutate(degree = igraph::degree(graph_network),
             group = dplyr::case_when(name %in% gene_symbol == TRUE ~ "Query Gene", 
                                      name %in% dep_network_top == TRUE ~ "Positive",
                                      name %in% dep_network_bottom == TRUE ~ "Negative"),
             group = as_factor(group), 
             group = fct_relevel(group, c("Query Gene", "Positive", "Negative")))  %>% #you don't end up with "connected" in a multi-gene list
      arrange(group) 
  }
  
  links <- graph_network %>%
    activate(edges) %>% # %E>%
    as_tibble()
  
  # determine the nodes that have at least the minimum degree
  nodes_filtered <- nodes %>%
    filter(degree >= deg) %>%  #input$degree
    as.data.frame
  
  # filter the edge list to contain only links to or from the nodes that have the minimum or more degree
  links_filtered <- links %>%
    filter(to %in% nodes_filtered$id & from %in% nodes_filtered$id) %>%
    as.data.frame
  
  # re-adjust the from and to values to reflect the new positions of nodes in the filtered nodes list
  links_filtered$from <- match(links_filtered$from, nodes_filtered$id) - 1
  links_filtered$to <- match(links_filtered$to, nodes_filtered$id) - 1
  
  #check to see if setting degree removed all links; if so, then throws error, so this fills a dummy links_filtered df to plot only nodes
  if(nrow(links_filtered) == 0) {links_filtered <- tibble("from" = 0, "to" = 0, "r2" = 1, "origin" = "pos")}
  
  #check to see if query gene is missing; if so, then adds a dummy so it shows up on graph, but disconnected
  if(sum(str_detect(nodes_filtered$group, "Query Gene")) == 0){
    dummy <- tibble("id" = max(nodes_filtered$id) + 1, "name" = gene_symbol, "degree" = 1, "group" = "Query Gene")
    nodes_filtered <- bind_rows(nodes_filtered, dummy)
  }
  # new visNetwork code from here down
  
  # get the approved_name of each gene from the gene_summary table - will be added to nodes tibble for tooltip
  nameTable <- tibble(name=character())
  for(gene in nodes_filtered$name){
    newVal <- gene_summary %>% 
      dplyr::filter(approved_symbol==gene) %>% 
      dplyr::pull(approved_name)
    if(length(newVal)==0){
      nameTable <- add_row(nameTable, name = "No gene summary information available")# handles cases where the gene is not in the gene summary table
    } else{
      nameTable <- add_row(nameTable, name=newVal)
    }
  }
  
  # add title information (tooltip that appears on hover) and 
  nodes_filtered <- nodes_filtered %>% 
    dplyr::mutate(title=paste0("<center><p>", nodes_filtered$name,"<br>",nameTable$name ,'<br><a target="_blank" href="https://datadrivenhypothesis.com/?show=gene&query_type=gene&symbol=',nodes_filtered$name,'">Gene Link</a></p>'),
                  label = nodes_filtered$name ) %>% 
    dplyr::mutate(id=0:(dim(nodes_filtered)[1]-1))
  
  # colors used within the network
  queryGeneColor <-"rgba(237, 165, 85, 0.8)"
  positiveColor <-"rgba(173, 103, 125, 0.8)"
  negativeColor <- "rgba(12, 35, 50, 0.8)"
  connectedColor <- "rgba(84, 64, 151, 0.8)"
  borderColor <- "rgba(255, 255, 255, 0.8)"
  edgeColor <- "rgba(84, 84, 84, 1)"
  
  #testing out variable spring lengths based upon the strength of association
  # links_filtered <- links_filtered %>%
  #   mutate(length = 100/(1-abs(r2)))
  
  # make node size as a function of the degree
  nodes_filtered <- nodes_filtered %>% 
    mutate(value = degree)
  
  # Setup filename for visExport
  if(length(gene_symbol)==1){
    exportName <- paste0(gene_symbol, "_networkGraph")
  } else{
    exportName <- ""
    for(geneName in gene_symbol){
      exportName <- paste0(geneName, "_", exportName)
    }
    exportName <- paste0("custom_",exportName, "_networkGraph")
  }
  
  
  # build the network visualization
  visNetwork(nodes = nodes_filtered,
             edges = links_filtered) %>% 
    visOptions(highlightNearest = c("enabled" = T, "hover" = T)) %>% 
    visGroups(groupname = "Query Gene", color = c("background" = queryGeneColor, 'border' =borderColor, 'highlight' = queryGeneColor, 'hover' = queryGeneColor ), shape='dot', borderWidth = 2) %>%
    visGroups(groupname = "Positive", color = c("background" = positiveColor, 'border' = borderColor, 'highlight' = positiveColor, 'hover' = positiveColor), shape='dot', borderWidth = 2) %>%
    visGroups(groupname = "Negative", color = c("background" = negativeColor, 'border' = borderColor, 'highlight' = negativeColor, 'hover' = negativeColor), shape='dot', borderWidth = 2) %>%
    visGroups(groupname = "Connected", color = c("background" = connectedColor, 'border' = borderColor, 'highlight' = connectedColor, 'hover' = connectedColor), shape='dot', borderWidth = 2) %>%
    visLegend(position = "right", width = .25, zoom = F) %>% 
    visEdges(color = edgeColor, smooth = F) %>% 
    visNodes(scaling = c("min" = 10, "max" =20)) %>% 
    visExport(name = exportName) %>% 
    visPhysics(barnesHut = c("damping" = 0.11), timestep = 0.25)  # reducing the timestep reduces the jitteriness of the graph and can help stabilize it
  # visConfigure(enabled=TRUE)
  # visSave(network, file = "networkTest2.html", background = "white")
}  

# Test Cases
# make_graph_visNetwork(gene_symbol = "SDHA")
# make_graph_visNetwork(gene_symbol = "CS",threshold = 18, deg = 2)

# make_graph_visNetwork(gene_symbol = "EXOC7",threshold = 20)

# make_graph_visNetwork(gene_symbol = c("SDHA", "SDHB"))
#make_graph_visNetwork(gene_symbol = c("GSS", "SST"))

make_graph_report <- function(toptable_data = master_top_table, bottomtable_data = master_bottom_table, gene_symbol, threshold = 10, deg = 2) {
  #get dep_network object
  dep_network <- setup_graph(toptable_data, bottomtable_data, gene_symbol, threshold)
  
  #make some objs for below
  if(length(gene_symbol) == 1){
    dep_top <- make_top_table(toptable_data, gene_symbol) %>% slice(1:threshold) #redundant with above, but need these objs
    dep_bottom <- make_bottom_table(bottomtable_data, gene_symbol) %>% slice(1:threshold)
  } else {
    dep_network_top <- dep_network %>% filter(origin == "pos") %>% pull(y)
    dep_network_bottom <- dep_network %>% filter(origin == "neg") %>% pull(y)
  }
  
  #make graph
  graph_network <- tidygraph::as_tbl_graph(dep_network)
  if(length(gene_symbol) == 1){
    nodes <-  as_tibble(graph_network) %>%
      rowid_to_column("id") %>%
      mutate(degree = igraph::degree(graph_network),
             group = case_when(name %in% gene_symbol == TRUE ~ "Query Gene", 
                               name %in% dep_top$Gene == TRUE ~ "Positive",
                               name %in% dep_bottom$Gene == TRUE ~ "Negative",
                               TRUE ~ "Connected"), 
             group = as_factor(group), 
             group = fct_relevel(group, c("Query Gene", "Positive", "Negative", "Connected")))  %>%
      arrange(group)
  } else {
    nodes <-  as_tibble(graph_network) %>%
      rowid_to_column("id") %>%
      mutate(degree = igraph::degree(graph_network),
             group = dplyr::case_when(name %in% gene_symbol == TRUE ~ "Query Gene", 
                                      name %in% dep_network_top == TRUE ~ "Positive",
                                      name %in% dep_network_bottom == TRUE ~ "Negative"),
             group = as_factor(group), 
             group = fct_relevel(group, c("Query Gene", "Positive", "Negative")))  %>% #you don't end up with "connected" in a multi-gene list
      arrange(group) 
  }
  
  links <- graph_network %>%
    activate(edges) %>% # %E>%
    as_tibble()
  
  # determine the nodes that have at least the minimum degree
  nodes_filtered <- nodes %>%
    filter(degree >= deg) %>%  #input$degree
    as.data.frame
  
  # filter the edge list to contain only links to or from the nodes that have the minimum or more degree
  links_filtered <- links %>%
    filter(to %in% nodes_filtered$id & from %in% nodes_filtered$id) %>%
    as.data.frame
  
  #readjust
  links_filtered$from <- match(links_filtered$from, nodes_filtered$id)
  links_filtered$to <- match(links_filtered$to, nodes_filtered$id)
  
  #check to see if setting degree removed all links; if so, then throws error, so this fills a dummy links_filtered df to plot only nodes
  if(nrow(links_filtered) == 0) {links_filtered <- tibble("from" = 1, "to" = 1, "r2" = 1, "origin" = "pos")}
  
  colors <- c("#EDA555","#AD677D", "#0C2332", "#544097")
  
  #check to see if query gene is missing; if so, then adds a dummy so it shows up on graph, but disconnected
  if(sum(str_detect(nodes_filtered$group, "Query Gene")) == 0){
    dummy <- tibble("id" = max(nodes_filtered$id) + 1, "name" = gene_symbol, "degree" = 1, "group" = "Query Gene")
    nodes_filtered <- bind_rows(nodes_filtered, dummy)
    #colors <- c("#AD677D", "#0C2332", "#544097", "#EDA555") #no need to reset colors; do that in 'breaks' below
  }
  
  graph_network_ggraph <- tidygraph::tbl_graph(nodes = nodes_filtered, edges = links_filtered)
  
  graph_network_ggraph %>%
    ggraph::ggraph(layout = "auto") +
    geom_edge_fan() + #edge_width = aes(abs(r2)), alpha = 0.3
    geom_node_point(aes(size = degree, color = group), alpha = 0.8) +
    geom_node_label(aes(filter = group != "Connected", label = name), repel = TRUE) +
    scale_colour_manual(values = colors, breaks = c("Query Gene", "Positive", "Negative", "Connected")) +
    theme_graph(base_family = 'Helvetica') +
    guides(size = "none", color = guide_legend(""))
  
  # # create visNetwork graph and save as a html
  # # get the approved_name of each gene from the gene_summary table - will be added to nodes tibble for tooltip
  # nameTable <- tibble(name=character())
  # for(gene in nodes_filtered$name){
  #   newVal <- gene_summary %>% 
  #     dplyr::filter(approved_symbol%in%nodes_filtered$name)  %>% # think I can delete this row...
  #     dplyr::filter(approved_symbol==gene) %>% 
  #     dplyr::pull(approved_name)
  #   if(length(newVal)==0){
  #     nameTable <- add_row(nameTable, name = "No gene summary information available")
  #   } else{
  #     nameTable <- add_row(nameTable, name=newVal)
  #   }
  # }
  # 
  # # add title information (tooltip that appears on hover) and 
  # nodes_filtered <- nodes_filtered %>% 
  #   dplyr::mutate(title=paste0("<center><p>", nodes_filtered$name,"<br>",nameTable$name ,'<br><a target="_blank" href="https://datadrivenhypothesis.com/?show=gene&query_type=gene&symbol=',nodes_filtered$name,'">Gene Link</a></p>'),
  #                 label = nodes_filtered$name ) %>% 
  #   dplyr::mutate(id=0:(dim(nodes_filtered)[1]-1))
  # 
  # # colors used within the network
  # queryGeneColor <-"rgba(237, 165, 85, 0.8)"
  # positiveColor <-"rgba(173, 103, 125, 0.8)"
  # negativeColor <- "rgba(12, 35, 50, 0.8)"
  # connectedColor <- "rgba(84, 64, 151, 0.8)"
  # borderColor <- "rgba(255, 255, 255, 0.8)"
  # edgeColor <- "rgba(84, 84, 84, 1)"
  # 
  # #testing out variable spring lengths based upon the strenght of association
  # # links_filtered <- links_filtered %>%
  # #   mutate(length = 100/(1-abs(r2)))
  # 
  # # node size as a function of the degree
  # nodes_filtered <- nodes_filtered %>% 
  #   mutate(value = degree)
  # 
  # # Setup filename for visExport
  # if(length(gene_symbol)==1){
  #   exportName <- paste0(gene_symbol, "_networkGraph")
  # } else{
  #   exportName <- ""
  #   for(geneName in gene_symbol){
  #     exportName <- paste0(geneName, "_", exportName)
  #   }
  #   exportName <- paste0("custom_",exportName, "_networkGraph")
  # }
  # 
  # 
  # # build the network visualization
  # visNetwork(nodes = nodes_filtered,
  #            edges = links_filtered,height = "100%", width = "100%") %>% 
  #   visOptions(highlightNearest = c("enabled" = T, "hover" = T)) %>% 
  #   visGroups(groupname = "Query Gene", color = c("background" = queryGeneColor, 'border' =borderColor, 'highlight' = queryGeneColor, 'hover' = queryGeneColor ), shape='dot', borderWidth = 2) %>%
  #   visGroups(groupname = "Positive", color = c("background" = positiveColor, 'border' = borderColor, 'highlight' = positiveColor, 'hover' = positiveColor), shape='dot', borderWidth = 2) %>%
  #   visGroups(groupname = "Negative", color = c("background" = negativeColor, 'border' = borderColor, 'highlight' = negativeColor, 'hover' = negativeColor), shape='dot', borderWidth = 2) %>%
  #   visGroups(groupname = "Connected", color = c("background" = connectedColor, 'border' = borderColor, 'highlight' = connectedColor, 'hover' = connectedColor), shape='dot', borderWidth = 2) %>%
  #   visLegend(position = "right", width = .25, zoom = F) %>% 
  #   visEdges(color = edgeColor, smooth = F) %>% 
  #   visNodes(scaling = c("min" = 10, "max" =20)) %>% 
  #   visExport(name = exportName) %>% 
  #   visPhysics(barnesHut = c("damping" = 0.11))
  # # visSave(network, file = "networkTest2.html", background = "white")
}

#figure legend
graph_title <- "Network Graph."
graph_legend <- "Each point represents a single gene taken from the top associated genes with the query gene. Genes with only one connection were removed."
graph_legend_list <- "Each point represents one of the queried genes, and then the top and bottom associated genes with it. Genes with only one connection were removed."

# make_graph_report(gene_symbol = "TP53")
#make_graph_report(gene_symbol = c("SDHA", "SDHB"))
#make_graph_report(gene_symbol = c("GSS", "SST"))

# exporting example... Yay! it can be output as a 2D or an html :)
# nodes <- data.frame(id = 1:3, group = c("B", "A", "B"))
# edges <- data.frame(from = c(1,2), to = c(2,3))
# 
# visNetwork(nodes, edges)  %>% 
#   visGroups(groupname = "A", color = "red")  %>% 
#   visGroups(groupname = "B", color = "lightblue")  %>% 
#   visLegend() %>% visExport()  
#   
#   visNetwork(nodes, edges) <!-- %>% -->
#   visGroups(groupname = "A", color = "red") <!-- %>% -->
#   visGroups(groupname = "B", color = "lightblue") <!-- %>% -->
#   visLegend() <!-- %>% visExport(type = "jpeg", name = "export-network",  -->
#                                    float = "left", label = "Save network", background = "purple", style= "") 
# 
# # saving as an html file example
# nodes <- data.frame(id = 1:3, group = c("B", "A", "B"))
# edges <- data.frame(from = c(1,2), to = c(2,3))
# 
# network <- visNetwork(nodes, edges)
# network
# 
# network 
# 
# # same as
# visSave(network, file = "network.html", background = "black")
