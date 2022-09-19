#-------------------------------------------------------------#
#                                                             #
#    GEMLI - Gene Expression Memory-based Lineage Inference   #
#                                                             #
#-------------------------------------------------------------#

# Set the working directory
setwd("~/Desktop/R_code/Memory_Almut/")

library(igraph)

# Load the GEMLI package
library(GEMLI)

# Load example data matrix (qc-filtered and normalized)
load('GEMLI_example_data_matrix.RData')

# Load the processed barcode information for testing of the predictions
load('GEMLI_example_barcode_information.RData')

#-------------------------------------------------------------#

# Example 1, small lineages, mESC

GEMLI_items = list()
GEMLI_items[['gene_expression']] = data_matrix
GEMLI_items[['barcodes']] = lineage_dict_bc

# Identify cell lineages through repeated iterative clustering (this may take 2-3min)
GEMLI_items = predict_lineages(GEMLI_items)
GEMLI_items$prediction[1:5,15:19]

# Test assigned lineages based on cell barcoding
GEMLI_items = test_lineages(GEMLI_items)
GEMLI_items$testing_results

# Visualize test results
GEMLI_items = test_lineages(GEMLI_items, plot_results=T)

# Visualize assigned lineages
visualize_as_network(GEMLI_items, cutoff=90)
visualize_as_network(GEMLI_items, cutoff=50)

# Visualize assigned lineages in comparison to a ground truth (e.g. from cell barcoding)
visualize_as_network(GEMLI_items, cutoff=90, ground_truth=T)
visualize_as_network(GEMLI_items, cutoff=50, ground_truth=T)

# Extract lineage information
GEMLI_items = prediction_to_lineage_information(GEMLI_items, cutoff=50)
GEMLI_items$predicted_lineage_table[1:10,]

# Suggest trimming for lineages that are too big
suggest_network_trimming_to_size(GEMLI_items, max_size=2, cutoff=50)

# Trim lineages that are too big
GEMLI_items_post_processed = trim_network_to_size(GEMLI_items, max_size=2, cutoff=50)

# Visualize trimmed network
visualize_as_network(GEMLI_items_post_processed, cutoff=50)
visualize_as_network(GEMLI_items_post_processed, cutoff=50, ground_truth=T)

#-------------------------------------------------------------#

# Call gene with gene expression memory based on predictions
GEMLI_items = markers_from_lineages(GEMLI_items)
GEMLI_items$lineage_markers[1:10,]

# rstudioapi::executeCommand('closeProject')

#-------------------------------------------------------------#

# Mind that the identification of cell lineages through repeated iterative clustering is a stochastic process. It therefore produces varying results each time it is run. To generate robust the results run at least 100 repetitions (default) or more.

#-------------------------------------------------------------#

#

#-------------------------------------------------------------#
#-------------------------------------------------------------#
#-------------------------------------------------------------#










#-------------------------------------------------------------#
#                       Work in progress                      #
#-------------------------------------------------------------#

check_lineage_assignment <- function(lineage_info)
{
  if (is.matrix(lineage_info)) {tmp = lineage_info[,"clone.ID"]; names(tmp) = lineage_info[,"cell.ID"]; lineage_info = tmp; rm(tmp)}
  table_to_plot = table(table(lineage_info))
  barplot(table_to_plot, names=names(table_to_plot), log="y", xlab="lineage size", ylab="lineage with size x", col="cornflowerblue")
}
check_lineage_assignment(lineage_prediction_info)

markers_by_variation <- function(data_matrix, family_dict, valid_fam_sizes=(2:7), use_median=T)
{
  family_dict_filt = family_dict[intersect(colnames(data_matrix), names(family_dict))]
  valid_family_dict = family_dict_filt[as.character(family_dict_filt) %in% names(table(family_dict_filt))[table(family_dict_filt) %in% valid_fam_sizes]]
  family_center = matrix(NA, ncol=length(unique(valid_family_dict)), nrow=nrow(data_matrix)); colnames(family_center) = unique(valid_family_dict); rownames(family_center) = rownames(data_matrix)
  for (family in as.character(unique(valid_family_dict)))
  {
    if (use_median){family_center[,family] = rowMedians(data_matrix[,names(valid_family_dict)[valid_family_dict==family]])}
    else {family_center[,family] = rowMeans(data_matrix[,names(valid_family_dict)[valid_family_dict==family]])}
  }
  x = rowMeans(family_center); y = cv_sq(family_center); x = log2(x); y = log2(y)
  filter = is.na(x) | is.na(y) | is.infinite(x) | is.infinite(y) | (x==0); x=x[!filter]; y=y[!filter]; loess_means = loess(y ~ x, span=0.75, control=loess.control(surface="direct"))
  filter = names(which(!filter))
  family_center_variation = loess_means$residuals
  # family_center_variation = log2(cv_sq(family_center[filter,])) - predict(loess_means, log2(rowMeans(family_center[filter,])))
  return(sort(family_center_variation[filter], decreasing=T))
}

markers_by_correlation <- function(data_matrix, family_dict, valid_fam_sizes=(2:7))
{
  genes = rownames(data_matrix)
  family_dict_filt = family_dict[intersect(colnames(data_matrix), names(family_dict))]
  valid_family_dict = family_dict_filt[as.character(family_dict_filt) %in% names(table(family_dict_filt))[table(family_dict_filt) %in% valid_fam_sizes]]
  real_family_matrix = outer(valid_family_dict, valid_family_dict, FUN='=='); diag(real_family_matrix) = F #; real_family_matrix[upper.tri(real_family_matrix)] = F
  tmp = which(real_family_matrix, arr.ind=T)
  family_members_a = colnames(real_family_matrix)[tmp[,1]]; family_members_b = colnames(real_family_matrix)[tmp[,2]]
  family_correlation <- function(data_matrix_row) {return(cor(data_matrix_row[family_members_a], data_matrix_row[family_members_b], method="spearman"))}
  family_corr = apply(data_matrix, 1, family_correlation)
  return(family_corr[order(family_corr, decreasing=T)])
}

marker_exons = markers_by_variation(data_matrix, lineage_prediction_info)
marker_exons = markers_by_correlation(data_matrix, lineage_prediction_info)

visualize_as_network <- function(clustering_results_repeated, lineage_dict_bc, cutoff=70, display_orphan=F)
{
  base_colors = rep(c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695','#40004b','#762a83','#9970ab','#c2a5cf','#e7d4e8','#d9f0d3','#a6dba0','#5aae61','#1b7837','#00441b','#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3','#e0e0e0','#bababa','#878787','#4d4d4d','#1a1a1a'), 100)
  network_edges = as.matrix(clustering_results_repeated)
  network_edges[network_edges<cutoff] <- 0
  if (!display_orphan) {network_edges = network_edges[rowSums(network_edges)!=0,colSums(network_edges)!=0]}
  vertex_color = base_colors[rank(as.numeric(lineage_dict_bc[rownames(network_edges)]))]; vertex_color[is.na(vertex_color)] = "lightgrey"
  network_edges = network_edges/max(network_edges)
  network_graph = igraph::graph.adjacency(network_edges, mode="undirected", weighted=T)
  network_graph = igraph::set.vertex.attribute(network_graph, "name", value=(1:ncol(network_edges)))
  edge_width = igraph::edge_attr(network_graph, "weight"); max(edge_width); min(edge_width)
  edge_width = edge_width - min(edge_width); edge_width = edge_width / max(edge_width)
  edge_width = edge_width * 9; max(edge_width); min(edge_width); edge_width = edge_width + 1; max(edge_width); min(edge_width) # alternative * 4.5 / + 0.5
  igraph::plot.igraph(network_graph, vertex.size=3, vertex.label.cex=0.5, vertex.color=vertex_color, vertex.label.color="black", vertex.label.family="Arial", edge.color="grey", rescale=TRUE, edge.width=edge_width, vertex.label.dist=0.5, vertex.label.degree=pi*1.5)
}

library(igraph)
prediction_to_lineage_information <- function(clustering_results_repeated, cutoff, output_as_dict=F)
{
  network = (clustering_results_repeated >= cutoff)
  network_edges = as.matrix(network)
  network_graph = igraph::graph.adjacency(network_edges, mode="undirected", weighted=NULL)
  family_dict = igraph::clusters(network_graph)$membership
  if (output_as_dict) {return(family_dict)}
  else {family_table = cbind(names(family_dict), as.vector(family_dict)); colnames(family_table) = c("cell.ID", "clone.ID"); return(family_table)}
}

#-------------------------------------------------------------#
