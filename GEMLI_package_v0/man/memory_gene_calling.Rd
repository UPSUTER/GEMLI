\name{memory_gene_calling}
\alias{memory_gene_calling}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
memory_gene_calling
}
\description{
This function identifies memory genes (or lineage markers) based on gene expression variability across lineages (either from predictions or from ground truth) from single-cell RNA-sequencing data.
}
\usage{
memory_gene_calling(GEMLI_items, valid_lineage_sizes=2:5, use_median=T, ground_truth=F)
}
\arguments{
  \item{
  GEMLI_items}{'GEMLI_items' is a list of GEMLI inputs and outputs. To run 'memory_gene_calling' GEMLI_items must contain a 'gene_expression' element. If it is run on bracodes 'GEMLI_items' needs to contain a 'barcodes' element. To run 'memory_gene_calling' on predictions 'GEMLI_items' needs to contain a 'predicted_lineage_table' elment. It is also possible to run 'memory_gene_calling' on a symmetric or asymmetric lineage type of specific cell types. In this case, the GEMLI_items must contain a 'cell_fate_analysis' element. The ‘gene_expression’ element is a quality controlled and normalised gene expression matrix where rownames are genes (features) and column names are cell IDs (cell barcodes). The 'predicted_lineage_table' is generated using the function prediction_to_lineage_information. The 'barcodes' element is a named vector of lineage ground truth (names=cell.ID, value=clone.ID). The 'cell_fate_analysis' element is generated using the 'extract_cell_fate_lineages' function.  
  }
  \item{
  valid_lineage_sizes}{'valid_lineage_sizes' specifies the range of lineage sizes to be included. Depending on the question to be investigated it can be beneficial to either restrict this to small lineages or large lineages respectively. Default is small lineage from 2 to 5 cells (2:5).
  }
  \item{
  use_median}{'use_median' specifies whether the median of lineages should be used rather than the mean. This makes the approach more robust to outliers. Default is 'true'/'T'.
  }
  \item{
  use_barcodes}{'ground_truth' specifies whether to call memory genes on ground truth instead of predictions. Default is 'false'/'F'.
  }
  \item{
  cell_fate}{'cell_fate', when present, specifies to call memory genes on specific symmetric or asymmetric lineages. It is a vector of the 'cell.fate' in the 'cell_fate_analysis' GEMLI_items element of the lineages to be used for memory gene calling. Cell fates start with sym or asym, followed by the cell type, separated by an underscore. 
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
'memory_gene_calling' yields a table of gene names or IDs of potential memory genes (row names), as well as their variability ('var') across lineages and a p-values ('p'). This table is stored un the GEMLI_items list as element 'memory_genes'.
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
