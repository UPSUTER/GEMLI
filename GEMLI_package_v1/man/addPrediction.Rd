\name{addPrediction}
\alias{addPrediction}
\alias{addPrediction,GEMLI-method}
\title{Add Prediction Method for GEMLI Objects}
\description{
  A generic function to add or update the prediction matrix in a \code{GEMLI} object.
}
\usage{
  \S4method{addPrediction}{GEMLI}(object, prediction)
}
\arguments{
  \item{object}{
    An object of class \code{GEMLI}.
  }
  \item{prediction}{
    A matrix containing prediction data.
  }
}
\details{
  This method updates the \code{prediction} slot of the \code{GEMLI} object. It checks that \code{prediction} is a matrix.
}
\value{
  Returns the updated \code{GEMLI} object.
}
\examples{
  gemli_obj <- GEMLI(gene_expression = matrix(rnorm(100), nrow = 10))
  gemli_obj <- addPrediction(gemli_obj, prediction = matrix(rnorm(20), nrow = 2))
}

