setClass(
  "GEMLI",
  slots = list(
    gene_expression = "dgCMatrix",     # Assuming gene_expression is a matrix
    barcodes = "character",         # Optional, vector of barcodes
    prediction = "matrix",      # Optional, data frame for prediction results
    testing_results = "matrix"  # Optional, data frame for testing results
  ),
  prototype = list(                 # Default values for optional slots
    barcodes = character(0),
    prediction = matrix(),
    testing_results = matrix()
  )
)

GEMLI <- function(gene_expression, barcodes=character(0), prediction=NULL, testing_results=NULL) {
    if (!is.matrix(gene_expression) & class(gene_expression)[1]!='Seurat' & class(gene_expression)[1]!='dgCMatrix' & class(gene_expression)[1]!='dgeMatrix') {
        stop("gene_expression must be a matrix, a Seurat object or a dgCMatrix")
    }
    if (!is.null(barcodes) && !is.character(barcodes)) {
        stop("barcodes must be a character vector")
    }
    if (!is.null(prediction) && !is.matrix(prediction)) {
        stop("prediction must be a data frame")
    }
    if (!is.null(testing_results) && !is.matrix(testing_results)) {
        stop("testing_results must be a data frame")
    }

    # Provide default empty data frames if NULL
    prediction <- if (is.null(prediction)) matrix() else prediction
    testing_results <- if (is.null(testing_results)) matrix() else testing_results

    if (class(gene_expression)[1]=='Seurat') {
        if (!requireNamespace("Seurat", quietly = TRUE)) {
            stop("Seurat is required for this function. Please install and load it before using cell_fate_DEG_calling")
        } else {
            gene_expression <- Matrix::Matrix(GetAssayData(object=gene_expression, layer="counts"))
        }
    }
    if (class(gene_expression)[1]=='dgeMatrix' | is.matrix(gene_expression)) {
        gene_expression <- Matrix::Matrix(gene_expression, sparse=T)
    }

    new("GEMLI", gene_expression = gene_expression, barcodes = barcodes, prediction = prediction, 
        testing_results = testing_results)
}

# Define the custom show method for the GEMLI class
setMethod("show", "GEMLI", function(object) {
  cat(paste0("GEMLI object with ", nrow(object@gene_expression), " genes and ", ncol(object@gene_expression), " cells\n\n"))

  if(length(object@barcodes) > 1) {
    cat("Barcodes: Provided\n")
  } else {
    cat("Barcodes: Not provided\n")
  }

  if(nrow(object@prediction) > 1) {
    cat("Prediction: Available\n")
  } else {
    cat("Prediction: Not available\n")
  }

  if(nrow(object@testing_results) > 1) {
    cat("Testing Results: Available\n")
  } else {
    cat("Testing Results: Not available\n")
  }
})

# Method to add or update barcodes
setGeneric("addBarcodes", function(object, barcodes) standardGeneric("addBarcodes"))
setMethod("addBarcodes", "GEMLI", function(object, barcodes) {
  if (!is.character(barcodes)) {
    stop("barcodes must be a character vector")
  }
  validObject(object@barcodes <- barcodes)
  return(object)
})

# Method to add or update prediction
setGeneric("addPrediction", function(object, prediction) standardGeneric("addPrediction"))
setMethod("addPrediction", "GEMLI", function(object, prediction) {
  if (!is.matrix(prediction)) {
    stop("prediction must be a matrix")
  }
  validObject(object@prediction <- prediction)
  return(object)
})

# Method to add or update testing results
setGeneric("addTestingResults", function(object, testing_results) standardGeneric("addTestingResults"))
setMethod("addTestingResults", "GEMLI", function(object, testing_results) {
  if (!is.matrix(testing_results)) {
    stop("testing results must be a matrix")
  }
  validObject(object@testing_results <- testing_results)
  return(object)
})
