calculate_correlations <- function(data, fast = FALSE)
{
    tmp = t(apply(data, 1, rank)); tmp = as.matrix(tmp); tmp = tmp - rowMeans(tmp); tmp = tmp / sqrt(rowSums(tmp^2))
    if(fast == TRUE){
        r = fastCor(t(tmp), nSplit = 10, upperTri = TRUE, optBLAS = TRUE, verbose=FALSE)}
    else{
        r = tcrossprod(tmp)
    }
    diag(r) <- 0
    return(r)
}
