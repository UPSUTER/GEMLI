potential_markers <- function(data_matrix) # check
{
  means = rowMeans(data_matrix)
  variation = (apply(data_matrix, 1, sd,  na.rm=T) / apply(data_matrix, 1, mean,  na.rm=T))**2
  filter = names(which(!(is.na(means) | is.na(variation) | is.infinite(means) | is.infinite(variation) | (means==0)))); x=means[filter]; y=variation[filter]
  linear_fit = lm(log2(y) ~ log2(x))
  variation_residuals = residuals(linear_fit)
  means = means[filter]
  mean_quantiles = quantile(means, seq(0.01,1,0.01))
  variation_quantiles = quantile(variation_residuals, seq(0.01,1,0.01))
  memory_genes = names(which((means>=mean_quantiles[98]) | (means>=mean_quantiles[90] & variation_residuals>=variation_quantiles[40]) | (means>=mean_quantiles[80] & variation_residuals>=variation_quantiles[80]) |(means>=mean_quantiles[60] & variation_residuals>=variation_quantiles[90])))
  return(memory_genes)
}
