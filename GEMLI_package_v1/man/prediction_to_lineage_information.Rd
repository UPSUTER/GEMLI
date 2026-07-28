\name{prediction_to_lineage_information}
\alias{prediction_to_lineage_information}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
prediction_to_lineage_information
}
\description{
This transforms a cell by cell matrix of confidence scores as created by the 'predict_lineages' function into a table (or vector, if specified) of predicted lineages.
}
\usage{
prediction_to_lineage_information(GEMLI_items, cutoff, output_as_dict = T)
}
\arguments{
  \item{
  GEMLI_items}{'GEMLI_items' is a list of GEMLI inputs and outputs. To run 'test_lineages' it should contain a prediction matrix named 'prediction' that is generated and added to the items list by the function 'predict lineages'.
  }
  \item{
  cutoff}{'cutoff' specifies the confidence score at which a cell pair is considered to be part of the same lineage. High values (e.g. 70-100) provide high precision but lower sensitivity. Low values (e.g. 30-60) provide higher sensitivity but lower precision. The default value is 50.
  }
  \item{
  output_as_dict}{'output_as_dict' AlmutNeeded }
}
\value{
The output is a matrix that lists cell IDs in the first column and the predicted lineage in the second column. This is added to the 'GEMLI_items' under the name 'predicted_lineage_table'. It also generates and adds the result as a vector that contains the predicted lineage as values and the cell IDs as names under the name 'predicted_lineages'.
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
