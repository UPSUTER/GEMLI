\name{GEMLI-class}
\alias{GEMLI-class}
\title{Class "GEMLI"}
\description{
  The \code{GEMLI} class represents a data structure for storing gene expression data, optional barcodes, prediction data, and testing results.
}
\details{
  The \code{GEMLI} class contains the following slots:
  \itemize{
    \item{\code{gene_expression}:} A \code{Matrix} object representing the gene expression data. It can be a matrix of counts or a sparse matrix.
    \item{\code{barcodes}:} A character vector of barcodes associated with the cells. It defaults to an empty character vector if not provided.
    \item{\code{prediction}:} A matrix containing prediction results. It defaults to an empty matrix if not provided.
    \item{\code{testing_results}:} A matrix containing testing results. It defaults to an empty matrix if not provided.
  }

  This class is designed to support the analysis and manipulation of gene expression data along with associated metadata.
}
\section{Slots}{
  \describe{
    \item{\code{gene_expression}}{A \code{Matrix} object.}
    \item{\code{barcodes}}{A character vector.}
    \item{\code{prediction}}{A matrix.}
    \item{\code{testing_results}}{A matrix.}
  }
}
\examples{
  # Example of creating a GEMLI object
  gemli_obj <- GEMLI(gene_expression = matrix(rnorm(100), nrow = 10))
  
  # Displaying the object
  show(gemli_obj)
}

