\name{trim_network_to_size}
\alias{trim_network_to_size}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
trim_network_to_size
}
\description{
This function trims predicted lineages that exceed a maximum lineage size. In lineages that exceed this size, the weakest links (least likely shared lineages) are trimmed until the desired maximum size is reached.
}
\usage{
trim_network_to_size(GEMLI_items, max_size = 4, cutoff = 70)
}
\arguments{
  \item{
  GEMLI_items}{'GEMLI_items' is a list of GEMLI inputs and outputs. To run 'suggest_network_trimming_to_size' it should contain a prediction matrix named 'prediction' that is generated and added to the items list by the function 'predict lineages'. It also needs to contain a vector of lineages (values are the lineages and names are the cell IDs) from predictions named 'predicted_lineages'. This vector can be added to the 'GEMLI_items' using the 'prediction_to_lineage_information' function.
  }
  \item{
  max_size}{'max_size' specifies maximum size of lineages.
  }
  \item{
  cutoff}{'cutoff' specifies the confidence score at which a cell pair is considered to be part of the same lineage. High values (e.g. 70-100) provide high precision but lower sensitivity. Low values (e.g. 30-60) provide higher sensitivity but lower precision. Default value is 70.
  }
}
\value{
This function return a 'GEMLI_items' list in which the prediction matrix ('prediction'), the prediction table ('predicted_lineage_table') and the prediction vector ('predicted_lineages') have been trimmed to the specified size.
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
