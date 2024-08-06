quantify_clusters_iterative = function(data_matrix, marker_genes, N=2, fast=FALSE)
{
  iterate = T; i = 2
  genes = intersect(marker_genes, rownames(data_matrix)[rowMeans(data_matrix)>0])
  data_matrix = data_matrix[genes,]
  corr_expr_raw = calculate_correlations(t(data_matrix), fast=FALSE); corr_expr = (1 - corr_expr_raw)/2
  cell_clusters = data.matrix(matrix(0, nrow=ncol(data_matrix), ncol=1)); rownames(cell_clusters) = colnames(data_matrix)
  cell_clusters[,1] = rep(1, ncol(data_matrix))
  while (iterate)
  {
    cell_clusters = cbind(cell_clusters, rep(0,nrow(cell_clusters)))
    for (cluster in setdiff(unique(cell_clusters[,(i-1)]),0))
    {
      cells_in_cluster = rownames(cell_clusters)[cell_clusters[,(i-1)]==cluster]
      if (length(cells_in_cluster) >= 4) # this line ends the sub clustering # min of desired cluster size
      {
        correlation = mean((corr_expr_raw[cells_in_cluster,cells_in_cluster])[lower.tri(corr_expr_raw[cells_in_cluster,cells_in_cluster], diag=F)])
        corr_expr_subset = corr_expr[cells_in_cluster,cells_in_cluster]
        clustering = cutree(hclust(as.dist(corr_expr_subset), method = "ward.D2", ), k=N)
        cell_clusters[names(clustering),i] = as.vector(clustering) + max(c(0, cell_clusters[,i]), na.rm=T)
      }
      else {cell_clusters[cells_in_cluster,i] = 0}
    }
    if (sum(cell_clusters[,i], na.rm=T)==0) {iterate = F}
    i = i+1
  }
  cell_clusters[cell_clusters==0] = NA
  return(cell_clusters)
}
