library(igraph)
prediction_to_lineage_information <- function(GEMLI_items, cutoff=50, output_as_dict=T) # check
{
  lineage_predictions_matrix = GEMLI_items[["prediction"]]
  network = (lineage_predictions_matrix >= cutoff)
  network_edges = as.matrix(network)
  network_graph = igraph::graph.adjacency(network_edges, mode="undirected", weighted=NULL)
  family_dict = igraph::clusters(network_graph)$membership
  GEMLI_items[["predicted_lineages"]] = family_dict
  family_table = cbind(names(family_dict), as.vector(family_dict)); colnames(family_table) = c("cell.ID", "clone.ID")
  GEMLI_items[["predicted_lineage_table"]] = family_table
  return(GEMLI_items)
}
