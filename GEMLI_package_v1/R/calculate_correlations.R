calculate_correlations <- function(data)
{
    tmp = t(apply(data, 1, rank)); tmp = as.matrix(tmp); tmp = tmp - rowMeans(tmp); tmp = tmp / sqrt(rowSums(tmp^2))
    r = tcrossprod(tmp)
    diag(r) <- 0
    return(r)
}
