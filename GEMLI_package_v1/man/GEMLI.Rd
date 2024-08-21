\name{GEMLI}
\alias{GEMLI}
\title{Create a GEMLI Object}
\description{
  Initializes a \code{GEMLI} object with gene expression data, optional barcodes, prediction, and testing results.
}
\usage{
  GEMLI(gene_expression, barcodes = character(0), prediction = NULL, testing_results = NULL)
}
\arguments{
  \item{gene_expression}{
    A matrix of gene expression data, a \code{Seurat} object, or a \code{dgCMatrix}. 
    If a \code{Seurat} object is provided, the function will extract the count matrix.
  }
  \item{barcodes}{
    A character vector of barcodes. Default is an empty character vector.
  }
  \item{prediction}{
    A matrix containing prediction data. Default is \code{NULL}, which will be replaced by an empty matrix.
  }
  \item{testing_results}{
    A matrix containing testing results. Default is \code{NULL}, which will be replaced by an empty matrix.
  }
}
\details{
  This function validates inputs for \code{gene_expression}, \code{barcodes}, \code{prediction}, and \code{testing_results}. It converts the gene expression data to a \code{Matrix} object if necessary, and ensures that \code{Seurat} is available if required.

  If \code{prediction} or \code{testing_results} are \code{NULL}, they are replaced with empty matrices.
}
\value{
  Returns an object of class \code{GEMLI} with the provided or default values.
}
\seealso{
  \code{\link{GEMLI-class}} for details on the \code{GEMLI} class.
}
\examples{
  # Example usage of the GEMLI function
  gemli_obj <- GEMLI(gene_expression = matrix(rnorm(100), nrow = 10))
}

