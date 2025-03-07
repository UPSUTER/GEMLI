\name{cluster_stability_plot}
\alias{cluster_stability_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
cluster_stability_plot
}
\description{
This function calculates a cluster stability index for lineage predictions allowing for different cluster sizes and generates a plot of cluster stability index vs cluster size. For multicellular structures, the cluster size at which the cluster stability index plateaus allows to estimate the maximal size of the multicellular structures present in the data. This is the cluster size till which lineage predictions will have a high precision, and the cluster size at which recovery will be maximal. 
}
\usage{
cluster_stability_plot(GEMLI_items)
}
\arguments{
  \item{
  GEMLI_items}{GEMLI_items is a list of GEMLI inputs and outputs. To run 'cluster_stability_plot' it should contain  a prediction matrix named 'prediction_multiple_sizes' that is generated and added to the items list by the function 'predict lineages_multiple_sizes'.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
'cluster_stability_plot' yields a plot of cluster stability index vs cluster size based on which the size of multicellular structures present in the single-cell RNA-sequencing dataset can be estimated.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Almut Eisele and Marcel Tarbier
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
