\name{DEG_volcano_plot}
\alias{DEG_volcano_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
DEG_volcano_plot
}
\description{
This function plots a simple volcano plot for the differential expressed genes (DEG) called using the function 'cell_fate_DEG_calling'.
}
\usage{
DEG_volcano_plot(GEMLI_items, name1, name2)
}
\arguments{
  \item{
  GEMLI_items}{GEMLI_items is a list of GEMLI inputs and outputs. '. To run 'DEG_volcano_plot' it should contain a 'DEG' element. The 'DEG' element is the output of the 'cell_fate_DEG_calling' function. It is a dataframe with the columns 'p_val', 'avg_log2FC', 'pct1', 'pct.2', 'p_val-adj'.
  }
  \item{
  name1}{'name1' is a character vector specifying the first population of cells analysed for DEG calling. It will appears in the title and legend of the volcano plot. It should correspond to the 'ident1' parameter of the 'cell.fate_DEG_calling' function used to generate the 'GEMLI_items' 'DEG' element.
  }
  \item{
  name2}{'name2' is a character vector specifying the second population of cells analysed for DEG calling. See 'name1' parameter.
  }
}
\value{
'DEG_volcano_plot' plots a volcano plot for the DEG called using the function 'cell_fate_DEG_calling'. 
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
