\name{quantify_clusters_iterative}
\alias{quantify_clusters_iterative}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
quantify_clusters_iterative
}
\description{
This function clusters the input matrix repeatedly until a desired cluster size is reached.
}
\usage{
quantify_clusters_iterative(data_matrix, marker_genes, N=2)
}
\arguments{
  \item{
  data_matrix}{'data_matrix' is a quality controlled and normalized gene expression matrix where rownames are genes (features) and column names are cell IDs (samples).
  }
  \item{
  marker_genes}{'marker_genes' is a vector of gene names or IDs of potential lineage marker genes. It is automatically created in the 'predict_lineages' function or can be computed manually using the 'potential_markers' function.
  }
  \item{
  N}{'N' describes in how many branches the data matrix is split in each clustering step. Higher numbers speed up the clustering but can negaatively impact the result of the prediction. It is highly recommended to keep the default value, 2.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The output is a matrix that indicates which cells (rows) cluster together in each iteration (colums). It is used in the 'predict_lineages' function.
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
