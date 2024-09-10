suggest_network_trimming_to_size <- function(GEMLI_items, max_size=4, cutoff=70, max_edge_width=5, display_orphan=F, include_labels=F, ground_truth=F, layout_style="fr")
{
    if (class(GEMLI_items)=='list') {
        barcodes = GEMLI_items[['barcodes']]
    } else if (class(GEMLI_items)=='GEMLI') {
        barcodes = GEMLI_items@barcodes
    } else {
        stop('Object GEMLI_items should be either of class list or GEMLI')
    }
  lineage_predictions_matrix_original = GEMLI_items[["prediction"]]
  predicted_lineages_original = GEMLI_items[["predicted_lineages"]]
  GEMLI_items_in_trimming = GEMLI_items
  predicted_lineages = prediction_to_lineage_information(GEMLI_items, cutoff, output_as_dict=T)$predicted_lineages
  if (ground_truth) {lineage_dict = barcodes} else {lineage_dict = prediction_to_lineage_information(GEMLI_items, cutoff, output_as_dict=T)$predicted_lineages}
  while (sum(table(GEMLI_items_in_trimming$predicted_lineages)>max_size)!=0)
  { predicted_lineages = GEMLI_items_in_trimming$predicted_lineages; 
    lineage_predictions_matrix = GEMLI_items_in_trimming$prediction
    oversized_families = names(table(predicted_lineages))[table(predicted_lineages)>max_size]
    for (oversized_family in oversized_families)
    {
      oversized_familiy_scores = lineage_predictions_matrix[names(which(predicted_lineages==oversized_family)), names(which(predicted_lineages==oversized_family))]
      weakest_link = min(oversized_familiy_scores[oversized_familiy_scores!=0])
      oversized_familiy_scores[oversized_familiy_scores==weakest_link] <- 0
      lineage_predictions_matrix[names(which(predicted_lineages==oversized_family)), names(which(predicted_lineages==oversized_family))] = oversized_familiy_scores
    }
    GEMLI_items_in_trimming$prediction = lineage_predictions_matrix
    GEMLI_items_in_trimming$predicted_lineages = prediction_to_lineage_information(GEMLI_items_in_trimming, cutoff, output_as_dict=T)$predicted_lineages
  }
  # visualization
  par(mar=c(0,0,3.5,0))
  layout(mat = matrix(c(1, 2), nrow = 2, ncol = 1), heights = c(3, 1))
  base_colors = rep(c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695','#40004b','#762a83','#9970ab','#c2a5cf','#e7d4e8','#d9f0d3','#a6dba0','#5aae61','#1b7837','#00441b','#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3','#e0e0e0','#bababa','#878787','#4d4d4d','#1a1a1a'), 100)
  network_edges = as.matrix(lineage_predictions_matrix_original)
  network_edges[lineage_predictions_matrix[rownames(network_edges),rownames(network_edges)]==0] <- ((-1) * network_edges[lineage_predictions_matrix[rownames(network_edges),rownames(network_edges)]==0])
  network_edges[abs(network_edges)<cutoff] <- 0
  if (!display_orphan) {network_edges = network_edges[rowSums(network_edges)!=0,colSums(network_edges)!=0]}
  
  vertex_color = base_colors[rank(as.numeric(lineage_dict[rownames(network_edges)]), na.last="keep", ties.method="min")]; vertex_color[is.na(vertex_color)] = "white"
    network_edges = network_edges/max(network_edges)
    network_graph = igraph::graph.adjacency(network_edges, mode="undirected", weighted=T)
    network_graph = igraph::set.vertex.attribute(network_graph, "name", value=(1:ncol(network_edges)))
    
    edge_color = igraph::edge_attr(network_graph, "weight") < 0; edge_color[edge_color==T] <- 'red'; edge_color[edge_color==F] <- 'grey'
      igraph::edge_attr(network_graph, "weight") = abs(igraph::edge_attr(network_graph, "weight"))
      edge_width = igraph::edge_attr(network_graph, "weight"); max(edge_width); min(edge_width)
      edge_width = edge_width - min(edge_width); edge_width = edge_width / max(edge_width)
      edge_width = edge_width * (max_edge_width*0.9); max(edge_width); min(edge_width); edge_width = edge_width + (max_edge_width*0.1); max(edge_width); min(edge_width)
      if (!include_labels) {network_graph = igraph::set.vertex.attribute(network_graph, "name", value=rep("",ncol(network_edges)))}
      if (layout_style=='fr') {igraph::plot.igraph(network_graph, vertex.size=3, vertex.label.cex=0.5, layout=igraph::layout.fruchterman.reingold(network_graph), vertex.color=vertex_color, vertex.label.color="black", vertex.label.family="Arial", edge.color=edge_color, rescale=TRUE, edge.width=edge_width, vertex.label.dist=0.5, vertex.label.degree=pi*1.5)}
      if (layout_style=='kk') {igraph::plot.igraph(network_graph, vertex.size=3, vertex.label.cex=0.5, layout=igraph::layout.kamada.kawai(network_graph), vertex.color=vertex_color, vertex.label.color="black", vertex.label.family="Arial", edge.color=edge_color, rescale=TRUE, edge.width=edge_width, vertex.label.dist=0.5, vertex.label.degree=pi*1.5)}
      # Title
      title(main = paste0("Prediction at confidence level ", cutoff, "\nsuggested trimming to size ", max_size, "\n"))
      # Edge_color_and_width
      plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
      legend("top", legend = c("Keep", "Trim"), col = c("grey", "red"), bty = "n", lty=1:1, lwd=c(max(edge_width), max(edge_width)), title = "Suggested", title.adj =0, horiz=F, xpd=TRUE, inset=c(0,0),ncol=1)
      legend("topleft", legend = c(max(lineage_predictions_matrix_original), min(lineage_predictions_matrix_original[which(lineage_predictions_matrix_original>cutoff)])), col = c("black", "black"), bty = "n", lty=1:1, lwd=c(max(edge_width),(min(edge_width)+ (max_edge_width*0.1))), title = "Confidence", title.adj =0, horiz=F, xpd=TRUE, inset=c(.15,0),ncol=1)
      # Vertex_color
      if (include_labels==T) {legend("topright", legend=c("Color - prediction","Number - cell ID"), bty = "n", title = "Vertex ", title.adj =0.5, horiz=F, xpd=TRUE, inset=c(.05, 0),ncol=1)} else
        {legend("topright", legend=c("Color - prediction"), bty = "n", title = "Vertex ", title.adj =0.5, horiz=F, xpd=TRUE, inset=c(.05, 0),ncol=1)}
}
