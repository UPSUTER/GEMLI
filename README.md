# Gene expression memory for lineage identification in scRNA-seq datasets

X is an R package to predict cell lineages (cells with a common ancestor) from single cell RNA datasets and to call genes with a high gene expression memory on the predicted cell lineages. It is described in publicationX.

## Installation
To install the package simply download the folder "Package_test". Open the project "LineAGED.Rproj" in RStudio. Then click in the RStudio menu on "Build" and then "Install and Restart". You can then close the window with the project and use the library in R.

## Getting started

### Load and explore example data

### Call potential lineage markers
Identify marker genes are identified based on gene expression mean and variation.
The `identify_markers` function takes a quality-controlled and normalized single-cell gene expression matrix (rows = genes/features, colums = cells/samples) as input. It outputs a vector of gene names or identifiers, depending on the input.
> `markers = identify_markers(data_matrix)`

### Perform lineage prediction

### Refine lineage prediction

### Examine lineage prediction

### Call real lineage markers

## Citation
If you use the package, please cite X

## TEST: Editing test by Marcel
<p align="center">
  <img width="500" height="500" src="https://github.com/UPSUTER/Memory/blob/main/Network_for_GitHub.png">
</p>
