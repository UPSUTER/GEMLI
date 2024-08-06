\name{calculate_correlations}
\alias{calculate_correlations}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
calculate_correlations
}
\description{
This function provides a fast way to calculate Spearmans ranked correlation or Pearsons correlation.
}
\usage{
calculate_correlations(data, fast=FALSE)
}
\arguments{
  \item{
  data}{'data_matrix' is a gene expression matrix where rownames are genes (features) and column names are cell IDs (samples).
  }
  \item{
  fast}{'fast' = FALSE will calculate the Spearman rank correlation, fast = TRUE will make use of the package HiClimR for the calculation of a Pearson correlation. The calculation of the Pearson correlation is faster than the calculation of Spearman rank correlation, however precision of lineage predictions will be slighlty reduced using the Pearson correlation. The default value is FALSE. 
}
}
\value{
The output is a cell by cell matrix with each value representing the rho value of Spearmans ranked correlation.
}
\references{
A.S. Eisele*, M. Tarbier*, A.A. Dormann, V. Pelechano, D.M. Suter | "Gene-expression memory-based prediction of cell lineages from scRNA-seq datasets" | Nature Communications 15, 2744 (2024). https://doi.org/10.1038/s41467-024-47158-y
}
\author{
Marcel Tarbier and Almut Eisele
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
