\name{suggest_network_trimming_to_size}
\alias{suggest_network_trimming_to_size}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
suggest_network_trimming_to_size
}
\description{
This function visualizes the lineage prediction as a network and highlights edges that could be removed based on a maximum lineage size cutoff. In lineages that exceed this size, the weakest links (least likely shared lineages) are suggested to be trimmed until the desired maximum size is reached. Trimming suggestions are highlighted in red.
}
\usage{
suggest_network_trimming_to_size(GEMLI_items) # max_size=4, cutoff=70, max_edge_width=5, display_orphan=F, include_labels=T, ground_truth=F)
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
  \item{
  max_edge_with}{'max_edge_width' specifies the maximum width of edges in the network visualization. All edge weights above the defined 'cutoff' will be scaled between 0.1*'max_edge_with' and 'max_edge_with'. Default value is 5.
  }
  \item{
  display_orphan}{'display_orphan' defines whether cells without connections should be displayed. This commonly leads the network plot getting less readable. Therefore the suggestion and the default are 'false'/'F'.
  }
  \item{
  include_labels}{'include_labels' defines whether nodes should be numbered. Cell IDs are not shown for readability.
  }
  \item{
  ground_truth}{If the 'GEMLI_items' list contains a 'barcodes' vector with orthogonal lineage information 'ground_truth' can be set 'true'/'T' to color cells according to their lineage. Default is 'false'/'F'.
  }
   \item{
  layout_style}{Depending on the number of cells, and the size of lineages, different layout styles can improve readability. In the suggest_network_trimming_to_size two layout algorithms can be chosen: Fruchterman-Reingold ("fr") and Kamada-Kawai ("kk"). Default is "fr".
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
This function has no output. It is for visualization only.
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
