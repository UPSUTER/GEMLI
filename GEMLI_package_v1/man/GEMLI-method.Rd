\name{show,GEMLI-method}
\alias{show,GEMLI-method}
\title{Show Method for GEMLI Objects}
\description{
  Displays a summary of the \code{GEMLI} object.
}
\usage{
  \S4method{show}{GEMLI}(object)
}
\arguments{
  \item{object}{
    An object of class \code{GEMLI}.
  }
}
\details{
  The \code{show} method provides a summary of the \code{GEMLI} object, including the number of genes and cells in the gene expression matrix, and the availability of barcodes, predictions, and testing results.
}
\value{
  Prints a summary to the console. No return value.
}
\examples{
  gemli_obj <- GEMLI(gene_expression = matrix(rnorm(100), nrow = 10))
  show(gemli_obj)
}

