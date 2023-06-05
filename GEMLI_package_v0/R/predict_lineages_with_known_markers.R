predict_lineages_with_known_markers <- function(GEMLI_items, repetitions=100, sample_size=(2/3), desired_cluster_size=c(2,3), N=2, fast=FALSE)
{
  norm_data = norm_data = GEMLI_items[['gene_expression']]
  marker_genes = GEMLI_items[['known_markers']]
  results = data.matrix(matrix(0, nrow=ncol(norm_data), ncol=ncol(norm_data))); rownames(results) = colnames(norm_data); colnames(results) = colnames(norm_data)
  for (i in seq(1,repetitions))
  {
    marker_genes_sample = sample(intersect(marker_genes, rownames(norm_data)), round(length(intersect(marker_genes, rownames(norm_data)))*sample_size,0))
    cell_clusters = quantify_clusters_iterative(norm_data, marker_genes_sample, N=2, fast=FALSE)
    cell_clusters_unique_name = cell_clusters; for (colname in 1:ncol(cell_clusters)){cell_clusters_unique_name[!is.na(cell_clusters_unique_name[,colname]),colname] = paste0(colname,'_',cell_clusters_unique_name[!is.na(cell_clusters_unique_name[,colname]),colname])}
    clustersize_dict = table(cell_clusters_unique_name)
    smallest_clusters = names(clustersize_dict)[clustersize_dict %in% desired_cluster_size]
    best_prediction = data.matrix(matrix(F, nrow=ncol(norm_data), ncol=ncol(norm_data))); rownames(best_prediction) = colnames(norm_data); colnames(best_prediction) = colnames(norm_data)
    for (cluster in smallest_clusters){cells_in_cluster = rownames(best_prediction)[rowSums(cell_clusters_unique_name==cluster, na.rm=T)>0]; best_prediction[cells_in_cluster,cells_in_cluster] <- T}
    diag(best_prediction) = F
    results = results + best_prediction
  }
  GEMLI_items[["prediction"]] = results
  return(GEMLI_items)
}
