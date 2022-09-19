trim_network_to_size <- function(GEMLI_items, max_size=4, cutoff=70)
{
  GEMLI_items_in_trimming = GEMLI_items
  while (sum(table(GEMLI_items_in_trimming$predicted_lineages)>max_size)!=0)
  {
    predicted_lineages = GEMLI_items_in_trimming$predicted_lineages; lineage_predictions_matrix = GEMLI_items_in_trimming$prediction
    oversized_families = names(table(predicted_lineages))[table(predicted_lineages)>max_size]
    for (oversized_family in oversized_families)
    {
      oversized_familiy_scores = lineage_predictions_matrix[names(which(predicted_lineages==oversized_family)), names(which(predicted_lineages==oversized_family))]
      weakest_link = min(oversized_familiy_scores[oversized_familiy_scores!=0])
      oversized_familiy_scores[oversized_familiy_scores==weakest_link] <- 0
      lineage_predictions_matrix[names(which(predicted_lineages==oversized_family)), names(which(predicted_lineages==oversized_family))] = oversized_familiy_scores
    }
    GEMLI_items_in_trimming$prediction = lineage_predictions_matrix
    GEMLI_items_in_trimming = prediction_to_lineage_information(GEMLI_items_in_trimming, cutoff)
  }
  return(GEMLI_items_in_trimming)
}
