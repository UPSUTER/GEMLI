predict_lineages <- function(GEMLI_items, repetitions=100, sample_size=(2/3), desired_cluster_size=c(2,3), fast=FALSE) # check
{
  if (class(GEMLI_items)=='list') {
      data_matrix = GEMLI_items[['gene_expression']]
  } else if (class(GEMLI_items)=='GEMLI') {
      data_matrix = GEMLI_items@gene_expression
  } else {
      stop('Object GEMLI_items should be either of class list or GEMLI')
  }
  marker_genes = potential_markers(data_matrix)
  results = data.matrix(matrix(0, nrow=ncol(data_matrix), ncol=ncol(data_matrix))); rownames(results) = colnames(data_matrix); colnames(results) = colnames(data_matrix)
  my_inters <- intersect(marker_genes, rownames(data_matrix))
  set.seed(42)
  myfun <- function(results) {
      marker_genes_sample = sample(my_inters, round(length(my_inters)*sample_size,0))
      cell_clusters = quantify_clusters_iterative(data_matrix, marker_genes_sample, N=2, fast)
      #
      cell_clusters_unique_name = cell_clusters
      for (colname in 1:ncol(cell_clusters)){
          cell_clusters_unique_name[!is.na(cell_clusters_unique_name[,colname]),colname] = paste0(colname,'_',cell_clusters_unique_name[!is.na(cell_clusters_unique_name[,colname]),colname])
      }
      #
      clustersize_dict = table(cell_clusters_unique_name)
      smallest_clusters = names(clustersize_dict)[clustersize_dict %in% desired_cluster_size]
      best_prediction = data.matrix(matrix(F, nrow=ncol(data_matrix), ncol=ncol(data_matrix)))
      rownames(best_prediction) = colnames(data_matrix)
      colnames(best_prediction) = colnames(data_matrix)
      for (cluster in smallest_clusters){
          cells_in_cluster = rownames(best_prediction)[rowSums(cell_clusters_unique_name==cluster, na.rm=T)>0]
          best_prediction[cells_in_cluster,cells_in_cluster] <- T
      }
      diag(best_prediction) = F
      results = results + best_prediction
      results
  }
  results <- Reduce(function(x, y) myfun(x), rep(list(results), repetitions))
  if (class(GEMLI_items)=='list') {
      GEMLI_items[["prediction"]] = results
  } else if (class(GEMLI_items)=='GEMLI') {
      GEMLI_items <- addPrediction(GEMLI_items, results)
  }
  return(GEMLI_items)
}
