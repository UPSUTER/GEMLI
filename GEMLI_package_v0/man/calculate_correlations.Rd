\name{calculate_correlations}
\alias{calculate_correlations}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
calculate_correlations
}
\description{
This function provides a fast way to calculate Spearman's ranked correlation.
}
\usage{
calculate_correlations(data_matrix)
}
\arguments{
  \item{
  data_matrix}{'data_matrix' is a gene expression matrix where rownames are genes (features) and column names are cell IDs (samples).
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The output is a cell by cell matrix with each value representing the rho value of Spearman's ranked correlation.
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
