\name{extract_cell_fate_lineages}
\alias{extract_cell_fate_lineages}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
extract_cell_fate_lineages
}
\description{
This function extracts symmetric (with all members in one considered cell type) and asymmetric cell lineages (with members in two or more of considered cell types). The function generates the input for the cell_fate_DEG_calling function.
}
\usage{
extract_cell_fate_lineages(GEMLI, selection, unique=FALSE, threshold)
}
\arguments{
  \item{
  GEMLI}{GEMLI is a list of GEMLI inputs and outputs. To run 'extract_cell_type_lineages' it should contain a 'predicted_lineage_table' as well as a 'cell_type' table. The 'predicted_lineage_table' is generated using the function 'prediction_to_lineage_information'. The 'cell_type' table is a data frame with column 'cell.ID' and celltype'.
  }
  \item{
  selection}{'selection' specifies the cell types to be considered for the extraction of symmetric and asymmetric lineages. It is a vector of minimal two characters specifying the cell types to be considered.
  }
  \item{
  unique}{'unique'=TRUE specifies that extracted lineages should contain only the cell types given in the 'selection' parameter. If 'unique'=FALSE, also other cell types, not considered in the lineage selection, can be present in the extraccted lineages. The default value is FALSE.
  }
  \item{
  threshold}{'threshold' specifies the minimal percentage of a given cell type asymmetric lineages should contain in order to be considered. It is a vector of percentages (numbers) which give the percentages for different cell types in the order in which they are given in the 'selection' parameter. Threshold values for all cell types have to be met for an asymmetric lineage to be kept.
  }
}
\value{
'extract_cell_fate_lineages' yields a data frame which is added to the 'GEMLI_items' under the name 'cell_fate_analysis'. The data frame contains the column 'cell.ID', 'clone.ID', 'cell.type' and 'cell.fate'. the 'cell.fate' column contains the lable 'asym' for selected asymmetric lineages and 'sym' for selected symmetric lineages, followed by the cell type of the specific cell, separated by an underscore (e.g. 'sym_DCIS', or 'asym_inv_tumor'). The function generates the input for the function 'cell_fate_DEG_calling'.
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
