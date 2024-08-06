\name{test_lineages}
\alias{test_lineages}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
test_lineages
}
\description{
This function tests the results of lineage assignments by comparing it to lineage assignments from cell barcoding.
}
\usage{
test_lineages(GEMLI_items, valid_fam_sizes=(1:5), max_interval=100, plot_results=F)
}
\arguments{
  \item{
  GEMLI_items}{'GEMLI_items' is a list of GEMLI inputs and outputs. To run 'test_lineages' it should contain a prediction matrix named 'prediction' that is generated and added to the items list by the function 'predict lineages'. It also needs to contain a ground truth names 'barcodes' provided as a named vector (values are the lineages and names are the cell IDs).
  }
  \item{
  valid_fam_sizes}{'valid_fam_sizes' specifies a resonable range for lineage sizes.
  }
  \item{
  max_interval}{'max_interval' is the number of repetitions used in the 'predict_lineages' function. The default value is 100.
  }
  \item{
  plot_results}{'plot_results' specifies whether the results of the test are visualized (plotted).
  }
}
\value{
The output is a table that lists the number of false positives (FP) and true positives (TP), the precision (TP/PP, where PP is the number of predicted positives which is the sum of TP and FP), and the sensitivity (TP/P, where P is the number of real positives which is the sum of TP and FN - false negatives).
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
}
