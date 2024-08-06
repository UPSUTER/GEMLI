\name{cell_type_composition_plot}
\alias{cell_type_composition_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
cell_type_composition_plot
}
\description{
This function generates simple plots of the cell type composition of predicted or ground truth lineages.
}
\usage{
cell_type_composition_plot(GEMLI_items, ground_truth=F, cell_type_colors=F, type, intersections=NULL)
}
\arguments{
  \item{
  GEMLI_items}{GEMLI_items is a list of GEMLI inputs and outputs. To run 'cell_type_composition_plot' it should contain a 'barcodes' (composition of ground truth) or 'predicted_lineage_table' (composition of predicted lineages) element. The 'predicted_lineage_table' is the output of the 'prediction_to_lineage_info' function. 
  }
  \item{
  ground_truth}{==T/TRUE indicates that the composition of ground truth lineages is analyzed. If 'ground_truth'==F, the compositionof predicted lineages is analyzed. Default is F.
  }
  \item{
  cell_type_colors}{'cell_type_colors'==T/TRUE specifies that custom colors for every cell type stored in GEMLI_items 'cell_type_colors' elment should be used. Default is F.
  }
  \item{
  type}{'type' specifies which of three plots is generated. Type can be 'bubble', 'upsetR', or 'plain'. type='plain' will output a simple table of the number of lineages for different cell type combinations. type='upsetR' will generate an upsetR plot showing the number of lineages for different cell type combinations. type='bubble' will generate a bubble plot of the cell type composition of individual lineages. This is especially meaningful when analyzing multicellular structures. 
  }
  \item{
  intersections}{'intersections' AlmutNeeded}
}
\value{
'cell_type_composition_plot' yields one of three possible plot types (see 'type') specifying the cell type composition of ground truth or predicted lineages.
}
\references{
A.S. Eisele*, M. Tarbier*, A.A. Dormann, V. Pelechano, D.M. Suter | "Gene-expression memory-based prediction of cell lineages from scRNA-seq datasets" | Nature Communications 15, 2744 (2024). https://doi.org/10.1038/s41467-024-47158-y
}
\author{
Almut Eisele and Marcel Tarbier
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
