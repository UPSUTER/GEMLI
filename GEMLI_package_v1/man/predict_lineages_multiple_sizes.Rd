\name{predict_lineages_multiple_sizes}
\alias{predict_lineages_multiple_sizes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
predict_lineages_multiple_sizes
}
\description{
This function predicts cell lineages from single-cell RNA-sequencing data with sizes ranging from a minimal to a maximal value. The prediction of lineages is performed as for the ‘predict_lineages’ function not only for a single desired lineage size, but independently for all cluster sizes in between a minimal and maximal value. For each prediction, the predicted lineage information is generated as for the ‘prediction_to_lineage_information function’. The function generates the input for the function ‘cluster_stability_plot’, which allows to estimate the size of multicellular structures present in the single-cell RNA-sequencing data.
  }
\usage{
predict_lineages_multiple_sizes(GEMLI_items, repetitions=10, sample_size=(2/3), minimal_maximal_cluster_size=c(2,50), cutoff=5)
}
\arguments{
  \item{
  GEMLI_items}{GEMLI_items is a list of GEMLI inputs and outputs. To run 'predict_lineages' it should contain a gene expression matrix named 'gene_expression'. This is a quality controlled and normalized gene expression matrix where rownames are genes (features) and column names are cell IDs (samples).
  }
  \item{
  repetitions}{'repetitions' specifies how many times the input matrix will be clustered using random subsamples of potential markers. A higher number of iterations leads to more robust results. 10 iterations is considered to be the minimum, but 100 iterations are strongly recommended (default value). Runtime is linear with regard to the number of iterations.
  }
  \item{
  sample_size}{'sample_size' is a value between 0 and 1 and specifies the fraction of potential markers that are used in each clustering. values between 0.5 and 0.67 (default value) are recommended.
  }
  \item{
  minimal_maximal_cluster_size}{'minimal_maximal_cluster_size' gives the minimal and maximal number of cells per cluster for which independent lineage predictions are run. The input is a vector of two values (minimal, maximal), e.g. c(2,50). The maximal value chosen should correspond to a value close or above the maximal expected size of multicellular structures present in the single-cell RNA-sequencing data. The default value is c(2,50).
  }
  \item{
  cutoff}{'cutoff' specifies the confidence score at which a cell pair is considered to be part of the same lineage. High values (e.g. 70-100) provide high precision but lower sensitivity. Low values (e.g. 30-60) provide higher sensitivity but lower precision.
  }
}
\value{
 'predict_lineages_multiple_sizes' yields a cell by lineage size matrix containing lineage IDs. This matrix is added to the 'GEMLI_items' under the name 'prediction_multiple_sizes'. The function generates the input for the function 'cluster_stability_plot', which allows to estimate the size of multicellular structures present in the single-cell RNA-sequencing data
}
\references{
A.S. Eisele*, M. Tarbier*, A.A. Dormann, V. Pelechano, D.M. Suter | "Gene-expression memory-based prediction of cell lineages from scRNA-seq datasets" | Nature Communications 15, 2744 (2024). https://doi.org/10.1038/s41467-024-47158-y
}
\author{
Almut Eisele and Marcel Tarbier
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
