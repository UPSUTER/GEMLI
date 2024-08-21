\name{addTestingResults}
\alias{addTestingResults}
\alias{addTestingResults,GEMLI-method}
\title{Add Testing Results Method for GEMLI Objects}
\description{
  A generic function to add or update the testing results matrix in a \code{GEMLI} object.
}
\usage{
  \S4method{addTestingResults}{GEMLI}(object, testing_results)
}
\arguments{
  \item{object}{
    An object of class \code{GEMLI}.
  }
  \item{testing_results}{
    A matrix containing testing results.
  }
}
\details{
  This method updates the \code{testing_results} slot of the \code{GEMLI} object. It checks that \code{testing_results} is a matrix.
}
\value{
  Returns the updated \code{GEMLI} object.
}
\examples{
  gemli_obj <- GEMLI(gene_expression = matrix(rnorm(100), nrow = 10))
  gemli_obj <- addTestingResults(gemli_obj, testing_results = matrix(rnorm(30), nrow = 3))
}

