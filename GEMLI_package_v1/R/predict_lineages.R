predict_lineages <- function(GEMLI_items, repetitions=100, sample_size=(2/3), desired_cluster_size=c(2,3),
                             verbose=T) {
    if (class(GEMLI_items)=='list') {
        data_matrix = Matrix::Matrix(GEMLI_items[['gene_expression']])
    } else if (class(GEMLI_items)=='GEMLI') {
        data_matrix = Matrix::Matrix(GEMLI_items@gene_expression)
    } else {
        stop('Object GEMLI_items should be either of class list or GEMLI')
    }

    cat('Define marker genes')
    marker_genes = potential_markers(data_matrix)

    cat('\nCompute predictions')
    #function to apply in parallel
    myfun <- function(idx, seed, verbose) {
        set.seed(seed)
        if (verbose) cat('  sample genes')
        marker_genes_sample = sample(marker_genes, round(length(marker_genes) * sample_size, 0))
        cell_clusters <- quantify_clusters_iterative(data_matrix, marker_genes_sample, N=2, verbose)
        clustersize_dict = table(paste0(
                                unlist(lapply(1:ncol(cell_clusters), function(i) rep(i, nrow(cell_clusters)))), 
                                '_', as.numeric(cell_clusters)
                                ))
        smallest_clusters = names(clustersize_dict)[clustersize_dict %in% desired_cluster_size]
        smallest_clusters = smallest_clusters[!grepl('_NA$', smallest_clusters)]
        best_prediction = matrix(FALSE, nrow=ncol(data_matrix), ncol=ncol(data_matrix))
        rownames(best_prediction) = colnames(data_matrix)
        colnames(best_prediction) = colnames(data_matrix)
        if (verbose) cat(', loop clusters')
        cells_in_cluster = lapply(smallest_clusters, function(cluster) {
            col <- as.numeric(strsplit(cluster, '_')[[1]][[1]])
            val <- as.numeric(strsplit(cluster, '_')[[1]][[2]])
            ans = na.omit(rownames(best_prediction)[cell_clusters[, col] == val])
            ans
        })
        for (i in 1:length(cells_in_cluster)) {
            best_prediction[cells_in_cluster[[i]], cells_in_cluster[[i]]] = T
        }
        diag(best_prediction) = F
        best_prediction
    }

    #compute
    set.seed(42) #for reproducibility
    results_list <- vector('list', repetitions)
    for (i in 1:repetitions) {
        if (verbose) cat(paste0('\n  Repetition', i, ': '))
        results_list[[i]] <- myfun(i, 42 +i, verbose)
    }

    #combine results
    cat('\nCombine results')
    results = matrix(0, nrow=ncol(data_matrix), ncol=ncol(data_matrix))
    rownames(results) = colnames(data_matrix)
    colnames(results) = colnames(data_matrix)
    for (mat in results_list) {
        results <- results + mat
    }

    #store results
    if (class(GEMLI_items) == 'list') {
        GEMLI_items[["prediction"]] = results
    } else if (class(GEMLI_items) == 'GEMLI') {
        GEMLI_items <- addPrediction(GEMLI_items, results)
    }

    return(GEMLI_items)
}
