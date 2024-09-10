#potential_markers <- function(data_matrix) # check
#{
  #means = rowMeans(data_matrix)
  #variation = (apply(data_matrix, 1, sd,  na.rm=T) / apply(data_matrix, 1, mean,  na.rm=T))**2
  #filter = names(which(!(is.na(means) | is.na(variation) | is.infinite(means) | is.infinite(variation) | (means==0)))); x=means[filter]; y=variation[filter]
  #linear_fit = lm(log2(y) ~ log2(x))
  #variation_residuals = residuals(linear_fit)
  #means = means[filter]
  #mean_quantiles = quantile(means, seq(0.01,1,0.01))
  #variation_quantiles = quantile(variation_residuals, seq(0.01,1,0.01))
  #memory_genes = names(which((means>=mean_quantiles[98]) | (means>=mean_quantiles[90] & variation_residuals>=variation_quantiles[40]) | (means>=mean_quantiles[80] & variation_residuals>=variation_quantiles[80]) |(means>=mean_quantiles[60] & variation_residuals>=variation_quantiles[90])))
  #return(memory_genes)
#}

potential_markers <- function(data_matrix) {
    # Calculate row means for a sparse matrix
    means <- Matrix::rowMeans(data_matrix)

    # Calculate standard deviation for each row in a sparse matrix
    row_sd_sparse <- function(sparse_matrix) {
        row_means <- Matrix::rowMeans(sparse_matrix)
        row_squares_means <- Matrix::rowMeans(sparse_matrix^2)
        sqrt(row_squares_means - row_means^2)
    }
    sds <- row_sd_sparse(data_matrix)

    # Calculate coefficient of variation squared
    variation <- (sds / means)^2

    # Filter out rows with NA, Inf, or means == 0
    filter <- which(!is.na(means) & !is.na(variation) & !is.infinite(means) & !is.infinite(variation) & means != 0)
    x <- means[filter]
    y <- variation[filter]

    # Perform linear regression on the log2-transformed data
    linear_fit <- lm(log2(y) ~ log2(x))

    # Get residuals from the linear model
    variation_residuals <- residuals(linear_fit)

    # Filter means and variation residuals
    means <- means[filter]

    # Calculate quantiles
    mean_quantiles <- quantile(means, seq(0.01, 1, 0.01))
    variation_quantiles <- quantile(variation_residuals, seq(0.01, 1, 0.01))

    # Identify memory genes based on quantiles
    memory_genes <- names(which(
                                (means >= mean_quantiles[98]) |
                                    (means >= mean_quantiles[90] & variation_residuals >= variation_quantiles[40]) |
                                    (means >= mean_quantiles[80] & variation_residuals >= variation_quantiles[80]) |
                                    (means >= mean_quantiles[60] & variation_residuals >= variation_quantiles[90])
                                ))

    return(memory_genes)
}

