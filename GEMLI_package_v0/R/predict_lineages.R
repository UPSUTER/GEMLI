predict_lineages <- function(GEMLI_items, repetitions=100, sample_size=(2/3), desired_cluster_size=c(2,3), N=2, fast=TRUE) # check
{
  data_matrix = GEMLI_items[['gene_expression']]
  marker_genes = potential_markers(data_matrix)
  results = data.matrix(matrix(0, nrow=ncol(data_matrix), ncol=ncol(data_matrix))); rownames(results) = colnames(data_matrix); colnames(results) = colnames(data_matrix)
  for (i in seq(1,repetitions))
  {
    marker_genes_sample = sample(intersect(marker_genes, rownames(data_matrix)), round(length(intersect(marker_genes, rownames(data_matrix)))*sample_size,0))
    cell_clusters = quantify_clusters_iterative(data_matrix, marker_genes_sample, N=2, fast)
    cell_clusters_unique_name = cell_clusters; for (colname in 1:ncol(cell_clusters)){cell_clusters_unique_name[!is.na(cell_clusters_unique_name[,colname]),colname] = paste0(colname,'_',cell_clusters_unique_name[!is.na(cell_clusters_unique_name[,colname]),colname])}
    clustersize_dict = table(cell_clusters_unique_name)
    smallest_clusters = names(clustersize_dict)[clustersize_dict %in% desired_cluster_size]
    best_prediction = data.matrix(matrix(F, nrow=ncol(data_matrix), ncol=ncol(data_matrix))); rownames(best_prediction) = colnames(data_matrix); colnames(best_prediction) = colnames(data_matrix)
    for (cluster in smallest_clusters){cells_in_cluster = rownames(best_prediction)[rowSums(cell_clusters_unique_name==cluster, na.rm=T)>0]; best_prediction[cells_in_cluster,cells_in_cluster] <- T}
    diag(best_prediction) = F
    results = results + best_prediction
  }
  GEMLI_items[["prediction"]] = results
  return(GEMLI_items)
}
