visualize_as_network <- function(GEMLI_items, cutoff=70, max_edge_with=5, display_orphan=F, include_labels=T, ground_truth=F, highlight_FPs=T)
{
  lineage_predictions_matrix = GEMLI_items[["prediction"]]
  if (ground_truth) {lineage_dict = GEMLI_items[['barcodes']]} else {lineage_dict = prediction_to_lineage_information(GEMLI_items, cutoff, output_as_dict=T)$predicted_lineages}
  base_colors = rep(c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695','#40004b','#762a83','#9970ab','#c2a5cf','#e7d4e8','#d9f0d3','#a6dba0','#5aae61','#1b7837','#00441b','#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3','#e0e0e0','#bababa','#878787','#4d4d4d','#1a1a1a'), 100)
  network_edges = as.matrix(lineage_predictions_matrix)
  network_edges[network_edges<cutoff] <- 0
  if (!display_orphan) {network_edges = network_edges[rowSums(network_edges)!=0,colSums(network_edges)!=0]}
  if ((ground_truth==T) & (highlight_FPs==T))
  {
    real_matrix = outer(lineage_dict[rownames(network_edges)], lineage_dict[rownames(network_edges)], "=="); real_matrix[is.na(real_matrix)] <- T
    rownames(real_matrix) = rownames(network_edges); colnames(real_matrix) = rownames(network_edges)
    network_edges[real_matrix[rownames(network_edges),rownames(network_edges)]==F] <- ((-1) * network_edges[real_matrix[rownames(network_edges),rownames(network_edges)]==F])
  }
  vertex_color = base_colors[rank(as.numeric(lineage_dict[rownames(network_edges)]), na.last="keep", ties.method="min")]; vertex_color[is.na(vertex_color)] = "white"
  network_edges = network_edges/max(network_edges)
  network_graph = igraph::graph.adjacency(network_edges, mode="undirected", weighted=T)
  network_graph = igraph::set.vertex.attribute(network_graph, "name", value=(1:ncol(network_edges)))
  edge_color = igraph::edge_attr(network_graph, "weight") < 0; edge_color[edge_color==T] <- 'red'; edge_color[edge_color==F] <- 'grey'
  igraph::edge_attr(network_graph, "weight") = abs(igraph::edge_attr(network_graph, "weight"))
  edge_width = igraph::edge_attr(network_graph, "weight"); max(edge_width); min(edge_width)
  edge_width = edge_width - min(edge_width); edge_width = edge_width / max(edge_width)
  edge_width = edge_width * (max_edge_with*0.9); max(edge_width); min(edge_width); edge_width = edge_width + (max_edge_with*0.1); max(edge_width); min(edge_width)
  if (!include_labels) {network_graph = igraph::set.vertex.attribute(network_graph, "name", value=rep("",ncol(network_edges)))}
  igraph::plot.igraph(network_graph, vertex.size=3, vertex.label.cex=0.5, vertex.color=vertex_color, vertex.label.color="black", vertex.label.family="Arial", edge.color=edge_color, layout=layout.fruchterman.reingold(network_graph), rescale=TRUE, edge.width=edge_width, vertex.label.dist=0.5, vertex.label.degree=pi*1.5)
}
