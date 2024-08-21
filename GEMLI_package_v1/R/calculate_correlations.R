calculate_correlations <- function(data)
{
    data <- t(as.matrix(data))
    tmp = t(apply(data, 1, rank))
    tmp = as.matrix(tmp)
    tmp = tmp - matrixStats::rowMeans2(tmp)
    tmp = tmp / sqrt(matrixStats::rowSums2(tmp^2))
    r = tcrossprod(tmp)
    diag(r) <- 0
    return(r)
}
