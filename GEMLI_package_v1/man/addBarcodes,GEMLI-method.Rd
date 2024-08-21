\name{addBarcodes,gemli-method}
\alias{addBarcodes,GEMLI-method}
\title{Add Barcodes Method for GEMLI Objects}
\description{
  Adds or updates the barcodes in a \code{GEMLI} object.
}
\usage{
  \S4method{addBarcodes}{GEMLI}(object, barcodes)
}
\arguments{
  \item{object}{
    An object of class \code{GEMLI}.
  }
  \item{barcodes}{
    A character vector of barcodes to be added or updated.
  }
}
\details{
  This method updates the \code{barcodes} slot of the \code{GEMLI} object. It checks that \code{barcodes} is a character vector.
}
\value{
  Returns the updated \code{GEMLI} object.
}
\examples{
  gemli_obj <- GEMLI(gene_expression = matrix(rnorm(100), nrow = 10))
  gemli_obj <- addBarcodes(gemli_obj, barcodes = c("barcode1", "barcode2"))
}

