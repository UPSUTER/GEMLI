visualize_as_network <- function(GEMLI_items, cutoff=70, max_edge_width=5, display_orphan=F, include_labels=F, ground_truth=F, highlight_FPs=F, layout_style='fr', cell_type_colors=F, title.cex=1, legend.cex=1)
{
    if (class(GEMLI_items)=='list') {
        barcodes = GEMLI_items[['barcodes']]
    } else if (class(GEMLI_items)=='GEMLI') {
        barcodes = GEMLI_items@barcodes
    } else {
        stop('Object GEMLI_items should be either of class list or GEMLI')
    }
    par(mar=c(1,1,5,1))
    if (cell_type_colors) {layout(mat = matrix(c(1, 2, 3, 0), nrow = 2, ncol = 2), heights = c(3, 1), widths = c(3,1))} else {layout(mat = matrix(c(1, 2), nrow = 2, ncol = 1), heights = c(3, 1))}
    # title
    lineage_predictions_matrix = GEMLI_items[["prediction"]]
    if (ground_truth) {lineage_dict = barcodes} else {lineage_dict = prediction_to_lineage_information(GEMLI_items, cutoff, output_as_dict=T)$predicted_lineages}
    base_colors = rep(c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695','#40004b','#762a83','#9970ab','#c2a5cf','#e7d4e8','#d9f0d3','#a6dba0','#5aae61','#1b7837','#00441b','#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3','#e0e0e0','#bababa','#878787','#4d4d4d','#1a1a1a'), 100)
    if (cell_type_colors) { if (length(GEMLI_items[['cell_type_color']])!=0){
        base_colors = GEMLI_items[['cell_type_color']]$color 
    } else {
        cell.type<-unique(GEMLI_items[['cell_type']]$cell.type)
        color<-base_colors[rank(unique(GEMLI_items[['cell_type']]$cell.type))]
        GEMLI_items[['cell_type_color']] = data.frame(cell.type, color)
    }} else {base_colors = base_colors
    }
    network_edges = as.matrix(lineage_predictions_matrix)
    network_edges[network_edges<cutoff] <- 0
    if (!display_orphan) {network_edges = network_edges[rowSums(network_edges)!=0,colSums(network_edges)!=0]}
    if ((ground_truth==T) & (highlight_FPs==T))
    { real_matrix = outer(lineage_dict[rownames(network_edges)], lineage_dict[rownames(network_edges)], "=="); real_matrix[is.na(real_matrix)] <- T
    rownames(real_matrix) = rownames(network_edges); colnames(real_matrix) = rownames(network_edges)
    network_edges[real_matrix[rownames(network_edges),rownames(network_edges)]==F] <- ((-1) * network_edges[real_matrix[rownames(network_edges),rownames(network_edges)]==F])
    }
    if (layout_style=='grid')
    { lineage_dict_filt = lineage_dict[rownames(network_edges)]
    family_order = names(sort(table(lineage_dict_filt), decreasing=T))
    cell_order = c(); for (family in family_order){cell_order = c(cell_order, names(lineage_dict_filt[lineage_dict_filt==family]))}
    network_edges = network_edges[cell_order, cell_order]
    }
    # Here is assignment of colors to vertex
    if (cell_type_colors) {vertex_color = GEMLI_items[['cell_type_color']]$color[match(GEMLI_items[['cell_type']]$cell.type[match(names(lineage_dict[rownames(network_edges)]), GEMLI_items[['cell_type']]$cell.ID)], GEMLI_items[['cell_type_color']]$cell.type)]} else
        {vertex_color = base_colors[rank(as.numeric(lineage_dict[rownames(network_edges)]), na.last="keep", ties.method="min")]}
    vertex_color[is.na(vertex_color)] = "white"
    network_edges = network_edges/max(network_edges)
    network_graph = igraph::graph.adjacency(network_edges, mode="undirected", weighted=T)
    network_graph = igraph::set.vertex.attribute(network_graph, "name", value=(1:ncol(network_edges)))
    edge_color = igraph::edge_attr(network_graph, "weight") < 0; edge_color[edge_color==T] <- 'red'; edge_color[edge_color==F] <- 'grey'
    igraph::edge_attr(network_graph, "weight") = abs(igraph::edge_attr(network_graph, "weight"))
    edge_width = igraph::edge_attr(network_graph, "weight");
    edge_width = edge_width - min(edge_width); edge_width = edge_width / max(edge_width) # get it to a 0-1 percentage value of the maximal present value
    edge_width = edge_width * (max_edge_width*0.9); edge_width = edge_width + (max_edge_width*0.1); 
    if (!include_labels) {network_graph = igraph::set.vertex.attribute(network_graph, "name", value=rep("",ncol(network_edges)))}
    if (layout_style=='fr') {plot.igraph(network_graph, vertex.size=3, vertex.label.cex=0.5, vertex.color=vertex_color, vertex.label.color="black", vertex.label.family="Arial", edge.color=edge_color, layout=igraph::layout.fruchterman.reingold(network_graph), rescale=TRUE, edge.width=edge_width, vertex.label.dist=0.5, vertex.label.degree=pi*1.5)}
    if (layout_style=='kk') {plot.igraph(network_graph, vertex.size=3, vertex.label.cex=0.5, vertex.color=vertex_color, vertex.label.color="black", vertex.label.family="Arial", edge.color=edge_color, layout=igraph::layout.kamada.kawai(network_graph), rescale=TRUE, edge.width=edge_width, vertex.label.dist=0.5, vertex.label.degree=pi*1.5)}
    if (layout_style=='grid') {plot.igraph(network_graph, vertex.size=3, vertex.label.cex=0.5, vertex.color=vertex_color, vertex.label.color="black", vertex.label.family="Arial", edge.color=edge_color, layout=igraph::layout.grid(network_graph), rescale=TRUE, edge.width=edge_width, vertex.label.dist=0.5, vertex.label.degree=pi*1.5)}
    # Title
    title(main = paste0("Prediction at confidence level ", cutoff, "\n "), cex.main=title.cex)
    # Edge_color_and_width
    plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
    if ((ground_truth==T) & (highlight_FPs==T)) {legend("bottom", legend = c("Correct", "False"), col = c("grey", "red"), bty = "n", lty=1:1, lwd=c(max(edge_width), max(edge_width)), title = "Prediction", title.adj =0, horiz=F, xpd=TRUE, ncol=1, cex=legend.cex)}
    legend("bottomleft", legend = c(max(lineage_predictions_matrix), min(lineage_predictions_matrix[which(lineage_predictions_matrix>cutoff)])), col = c("black", "black"), bty = "n", lty=1:1, lwd=c(max(edge_width),(min(edge_width)+ (max_edge_width*0.1))), title = "Confidence", title.adj =0, horiz=F, xpd=TRUE, ncol=1, cex=legend.cex)
    # Vertex_color
    if ((ground_truth==T) & (cell_type_colors==F) & (include_labels==F)) {legend("bottomright", legend=c("Color - ground truth","White - no ground truth"), bty = "n", title = "Color ", title.adj =0.5, horiz=F, xpd=TRUE, inset=c(.05, 0),ncol=1, cex=legend.cex)}
    if ((ground_truth==F) & (cell_type_colors==F) & (include_labels==F)) {legend("bottom", legend=c("prediction",""), bty = "n", title = "  Color by", title.adj =0.5, horiz=F, xpd=TRUE,ncol=1, cex=legend.cex)}
    if ((ground_truth==T) & (cell_type_colors==F) & (include_labels==T)) {legend("topright", legend=c("Color - ground truth","White - no ground truth","Number - cell ID"), bty = "n", title = "Vertex ", title.adj =0.5, horiz=F, xpd=TRUE, inset=c(.05, 0),ncol=1, cex=legend.cex)}
    if ((ground_truth==F) & (cell_type_colors==F) & (include_labels==T)) {legend("top", legend=c("Color - prediction","Number - cell ID"), bty = "n", title = "  Vertex", title.adj =0.5, horiz=F, xpd=TRUE, inset=c(0, 0),ncol=1, cex=legend.cex)}
    if (cell_type_colors) {
        plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
        legend("left", legend = GEMLI_items[['cell_type_color']]$cell.type, pch = 16, col = GEMLI_items[['cell_type_color']]$color, title = "Cell type", bty = "o", horiz=F, xpd=TRUE, inset=c(0, 0),ncol=1, cex=legend.cex)}
    par(mar=c(5.1, 4.1, 4.1, 2.1))
}
