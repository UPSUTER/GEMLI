quantify_clusters_iterative <- function(data_matrix, marker_genes, N = 2, verbose=T) {
    iterate <- TRUE
    i <- 2
    if (verbose) cat(', get correlation')
    corr_expr_raw <- (1 - calculate_correlations(data_matrix[marker_genes, , drop=F]) ) / 2
    cell_clusters <- Matrix::Matrix(0, nrow = ncol(data_matrix), ncol = 1, sparse = F)
    rownames(cell_clusters) <- colnames(data_matrix)
    cell_clusters[, 1] <- rep(1, ncol(data_matrix))

    if (verbose) cat(', iterate')
    while (iterate) {
        cell_clusters <- cbind(cell_clusters, Matrix::Matrix(0, nrow = nrow(cell_clusters), sparse = F))
        unique_clusters <- setdiff(unique(cell_clusters[, (i - 1)]), 0)

        for (cluster in unique_clusters) {
            cluster_indices <- which(cell_clusters[, (i - 1)] == cluster)
            cells_in_cluster <- rownames(cell_clusters)[cluster_indices]

            if (length(cells_in_cluster) >= 4) {
                clustering <- cutree(fastcluster::hclust(as.dist(corr_expr_raw[cells_in_cluster, cells_in_cluster]), method = "ward.D2"), k = N)
                cell_clusters[cluster_indices, i] <- as.vector(clustering) + max(c(0, cell_clusters[, i]), na.rm = TRUE)
            }
        }
        dim(cell_clusters)

        if (sum(cell_clusters[, i], na.rm = TRUE) == 0) {
            iterate <- FALSE
        }
        i <- i + 1
    }
    lapply(c("corr_expr_raw", "unique_clusters", "cluster_indices", "cells_in_cluster", "clustering"), function(obj) if (exists(obj)) rm(list = obj))
    cell_clusters[cell_clusters == 0] <- NA
    return(cell_clusters)
}
