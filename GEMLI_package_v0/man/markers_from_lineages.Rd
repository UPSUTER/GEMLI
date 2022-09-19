\name{markers_from_lineages}
\alias{markers_from_lineages}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
markers_from_lineages
}
\description{
This function identifies lineage markers based on gene expression variability across lineages (either from predictions or from barcodes) from single-cell RNA-sequencing data.
}
\usage{
markers_from_lineages(GEMLI_items)
}
\arguments{
  \item{
  GEMLI_items}{'GEMLI_items' is a list of GEMLI inputs and outputs. To run 'markers_from_lineages' it should contain a prediction matrix named 'prediction' that is generated and added to the items list by the function 'predict lineages'. It also needs to contain a vector of lineages (values are the lineages and names are the cell IDs) either from predictions or barcodes. If it is run on predictions 'GEMLI_items' needs to contain a vector named 'predicted_lineages'. If it is run on bracodes 'GEMLI_items' needs to contain a vector named 'barcodes'.
  }
  \item{
  valid_lineage_sizes}{'valid_lineage_sizes' specifies the range of lineage sizes to be included. Depending on the question to be investigated it can be beenficial to either restrict this to small lineages or large lineages respectively. Default is small lineage from 2 to 5 cells (2:5).
  }
  \item{
  use_median}{'use_median' specifies whether the median of lineages should be used rather than the mean. This makes the approach more robust to outliers. Default is 'true'/'T'.
  }
  \item{
  use_barcodes}{'use_barcodes' specifies whether to call lineage markers on barcodes instead of predictions. Default is 'false'/'F'.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
'markers_from_lineages' yields a table of gene names or IDs of potential lineage marker genes (row names), as well as their variability ('var') across lineages and p-values ('p').
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Marcel Tarbier and Almut Eisele
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
