predict_lineages_multiple_sizes <- function(GEMLI_items, repetitions=10, sample_size=(2/3), minimal_maximal_cluster_size=c(2,50), N=2, cutoff=5) # check
{
  # split out the desired_cluster_size into vectors encompassing all combis between the min and max value
  # for each of the vectors then generate the lineage prediction and combine them into one table
  desired_sizes<-list()
  for (m in 2:minimal_maximal_cluster_size[2]){desired_sizes[[m-1]]<-c(1:m)}
  GEMLI_prediction_list<-list()
  for (j in 1:length(desired_sizes)) {
  progress((100*j)/length(desired_sizes))
  sub_desired_cluster_size<-desired_sizes[[j]]
  data_matrix = GEMLI_items[['gene_expression']]
  marker_genes = potential_markers(data_matrix)
  results = data.matrix(matrix(0, nrow=ncol(data_matrix), ncol=ncol(data_matrix))); rownames(results) = colnames(data_matrix); colnames(results) = colnames(data_matrix)
  for (i in seq(1,repetitions))
  {
    marker_genes_sample = sample(intersect(marker_genes, rownames(data_matrix)), round(length(intersect(marker_genes, rownames(data_matrix)))*sample_size,0))
    cell_clusters = quantify_clusters_iterative(data_matrix, marker_genes_sample, N=2)
    cell_clusters_unique_name = cell_clusters; for (colname in 1:ncol(cell_clusters)){cell_clusters_unique_name[!is.na(cell_clusters_unique_name[,colname]),colname] = paste0(colname,'_',cell_clusters_unique_name[!is.na(cell_clusters_unique_name[,colname]),colname])}
    clustersize_dict = table(cell_clusters_unique_name)
    
    # This is where the desired_cluster_size comes into play
    smallest_clusters = names(clustersize_dict)[clustersize_dict %in% sub_desired_cluster_size]
    best_prediction = data.matrix(matrix(F, nrow=ncol(data_matrix), ncol=ncol(data_matrix))); rownames(best_prediction) = colnames(data_matrix); colnames(best_prediction) = colnames(data_matrix)
    for (cluster in smallest_clusters){cells_in_cluster = rownames(best_prediction)[rowSums(cell_clusters_unique_name==cluster, na.rm=T)>0]; best_prediction[cells_in_cluster,cells_in_cluster] <- T}
    diag(best_prediction) = F
    results = results + best_prediction
  }
  # Transform the prediction matrix into a lineage info
  network = (results >= cutoff)
  network_edges = as.matrix(network)
  network_graph = igraph::graph.adjacency(network_edges, mode="undirected", weighted=NULL)
  prediction_dict = igraph::clusters(network_graph)$membership
  family_table = cbind(names(prediction_dict), as.vector(prediction_dict))
  colnames(family_table) = c("cell.ID", "clone.ID")
  GEMLI_prediction_list[[j]] = family_table
  }
  
  # combine in one dataframe that can be used for cluster stability calculation
  for(u in 1:length(GEMLI_prediction_list)){colnames(GEMLI_prediction_list[[u]]) <- c("cell.ID",paste0("K", u)) }
  GEMLI_prediction_multiple_sizes_output <- Reduce(function(x,y)merge(x,y,by="cell.ID"), GEMLI_prediction_list)
  GEMLI_items[['prediction_multiple_sizes']]<-GEMLI_prediction_multiple_sizes_output
  return(GEMLI_items)
}
