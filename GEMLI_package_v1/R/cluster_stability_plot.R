cluster_stability_plot <- function(GEMLI_items) # check
{
    data_matrix<-GEMLI_items[['prediction_multiple_sizes']]
    clustree<-clustree(data_matrix, prefix = "K")
    clustree<-clustree[["data"]][which(clustree[["data"]]$size != 1),]
    plot(clustree$K, clustree$sc3_stability, xlab="lineage size", ylab="clustree_stability_index")
}
