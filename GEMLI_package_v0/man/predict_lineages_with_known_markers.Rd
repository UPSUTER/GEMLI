\name{predict_lineages_with_known_markers}
\alias{predict_lineages_with_known_markers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
predict_lineages_with_known_markers
}
\description{
This function predicts cell lineages from from single-cell RNA-sequencing data when lineage markers are already knwon, e.g. from a barcoding experiment. Subsets of these genes are used to cluster cells until the desired cluster size is reached. This clustering is repeated many times for random subsets. The result is a cell by cell matrix that lists how many times each cell pair clustered together which translates into a confidence score. Cell pairs with high confidence scores are likely to be members of the same lineage.
}
\usage{
predict_lineages_with_known_markers(GEMLI_items, repetitions=100, sample_size=(2/3), desired_cluster_size=c(2,3), fast=FALSE)
}
\arguments{
  \item{
  GEMLI_items}{GEMLI_items is a list of GEMLI inputs and outputs. To run 'predict_lineages' it should contain a gene expression matrix named 'gene_expression'. This is a quality controlled and normalized gene expression matrix where rownames are genes (features) and column names are cell IDs (samples). In addition it needs to contain a vector of known marker genes names named 'known_markers'.
  }
  \item{
  repetitions}{'repetitions' specifies how many times the input matrix will be clustered using random subsamples of potential markers. A higher number of iterations leads to more robust results. 10 iterations is considered to be the minimum, but 100 iterations are strongly recommended (default value). Runtime is linear with regard to the number of iterations.
  }
  \item{
  sample_size}{'sample_size' is a value between 0 and 1 and specifies the fraction of potential markers that are used in each clustering. values between 0.5 and 0.67 (default value) are recommended.
  }
  \item{
  desired_cluster_size}{'desired_cluster_size' specifies the number of cells per cluster to be achieved in each clustering. The input is a lsit of values, e.g. c(2,3,4) or (2:4). The desired_cluster_size parameter should generally be small. Values between 2 and 4 are recommended (default).
  }
  \item{
  fast}{'fast' =TRUE uses the HiClimR package for calculating correlations. This makes predictions faster but reduces precision. The default value is FALSE.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
'predict_lineages' yields a cell by cell matrix containing confidence scores. Cell pairs with high confidence scores are likely to be members of the same lineage. This matrix is added to the 'GEMLI_items' under the name 'prediction'.
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
