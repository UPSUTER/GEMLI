\name{potential_marker}
\alias{potential_marker}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
potential_marker
}
\description{
This function identifies potential lineage markers based on mean gene expression and gene expression variability from single-cell RNA-sequencing data. It is part of the 'predict_lineages' function, but can also be called idependently on a quality controlled and normalized gene expression matrix where rownames are genes (features) and column names are cell IDs (samples).
}
\usage{
potential_marker(data_matrix)
}
\arguments{
  \item{
  data_matrix}{'data_matrix' is a quality controlled and normalized gene expression matrix where rownames are genes (features) and column names are cell IDs (samples).
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
'potential_marker' yields a vector of gene names or IDs of potential lineage marker genes. While genes are selected purely based on gene expression mean and variability, it has been shown that this approach enriches for genes with lineage specific gene expression profiles.
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
