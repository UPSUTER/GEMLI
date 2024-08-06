test_lineages <- function(GEMLI_items, valid_fam_sizes=(1:5), max_interval=100, plot_results=F)
{
  lineage_predictions_matrix = GEMLI_items[['prediction']]
  lineage_dict_bc = GEMLI_items[['barcodes']]
  valid_family_dict = lineage_dict_bc[as.character(lineage_dict_bc) %in% names(table(lineage_dict_bc))[table(lineage_dict_bc) %in% valid_fam_sizes]]
  cell_with_annotation = intersect(rownames(lineage_predictions_matrix), names(valid_family_dict))
  family_dict_filt = valid_family_dict[cell_with_annotation]
  real_family_matrix = outer(family_dict_filt[cell_with_annotation], family_dict_filt[cell_with_annotation], FUN='=='); diag(real_family_matrix) = F
  results_repeated_annotated = lineage_predictions_matrix[cell_with_annotation, cell_with_annotation]
  if (is.na(max_interval)) {intervals = unique(round(seq(0,1,0.1)*max(results_repeated_annotated),0))} else {intervals = unique(round(seq(0,1,0.1)*max_interval,0))}
  output_matrix = matrix(NA, ncol=4, nrow=length(intervals)); colnames(output_matrix) = c('precision','TP','FP','sensitivity'); rownames(output_matrix) = intervals
  for (interval in intervals){output_matrix[as.character(interval),1:3] = c(sum(real_family_matrix & (results_repeated_annotated>=interval)) / sum(results_repeated_annotated>=interval), sum(real_family_matrix & (results_repeated_annotated>=interval)), sum((!real_family_matrix) & (results_repeated_annotated>=interval)))}
  # replace this with an apply function
  output_matrix[,"sensitivity"] = output_matrix[,"TP"]/output_matrix["0","TP"]
  output_matrix = output_matrix[,c('TP','FP','precision','sensitivity')]
  if (plot_results)
  {
    par(mar=c(4.5, 4.5, 3.5, 4.5)); plot(as.numeric(rownames(output_matrix)), output_matrix[,"precision"], type="o", pch=16, lwd=3, col="darkred", xlab="confidence level", ylab="precision (red)", log="", main="testing lineage prediction", ylim=c(0,1)); par(new=T); plot(as.numeric(rownames(output_matrix)), output_matrix[,"sensitivity"], type="o", pch=16, lwd=3, axes=F, bty="n", xlab="", ylab="", col="grey", log="", ylim=c(0,1)); axis(side=4, at=pretty(range(output_matrix[,"sensitivity"]))); mtext("sensitivity (grey)", side=4, line=3)
  }
  GEMLI_items[['testing_results']] = output_matrix
  return(GEMLI_items)
}
