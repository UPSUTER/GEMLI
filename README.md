# Gene expression memory for lineage identification in scRNA-seq datasets

X is an R package to predict cell lineages (cells with a common ancestor) from single cell RNA datasets and to call genes with a high gene expression memory on the predicted cell lineages. It is described in publicationX.

## Installation
To install the package simply download the folder "Package_test". Open the project "LineAGED.Rproj" in RStudio. Then click in the RStudio menu on "Build" and then "Install and Restart". You can then close the window with the project and use the library in R.

## Getting started

### Call potential lineage markers
Identify marker genes are identified based on gene expression mean and variation.
The `identify_markers` function takes a quality-controlled and normalized single-cell gene expression matrix (rows = genes/features, colums = cells/samples) as input. It outputs a vector of gene names or identifiers, depending on the input.

```
> markers = identify_markers(data_matrix)
> markers[1:3]
[1] "ENSMUSG00000025915" "ENSMUSG00000048960" "ENSMUSG00000043716"
```

### Perform lineage prediction
Identify cell lineages through repeated iterative clustering (this may take 2-3min).
The `predict_lineages` function takes a quality-controlled and normalized single-cell gene expression matrix (rows = genes/features, colums = cells/samples) as input. It outputs a matrix of all cells against all cells with values corresponding to a confidence score that they are part of the same lineage. 

```
> lineage_predictions_matrix = predict_lineages(data_matrix)
> lineage_predictions_matrix[1:5,15:19]
                   AGAGAATAGGTCATAA-1 AGAGCAGCAAGTGATA-1 AGATGCTTCAAAGACA-1 AGGATCTGTATCGTTG-1 AGGGAGTAGACGATAT-1
AAACGAACAGGTGTGA-1                  0                  0                  0                  0                 14
AAAGGTAGTTGCTTGA-1                 87                  0                  0                  0                  0
AACCACAAGTTTGTCG-1                  0                  0                  0                  0                 39
AAGCCATGTTCCACGG-1                  0                  0                  0                  0                  0
AAGCGAGGTACGGCAA-1                  0                  0                  0                  0                  0
```

### Test lineage prediction
Test predicted lineages against lineages from cell barcoding.
The `test_lineage_prediction` function takes the output of `predict_lineages` as well as a the file 'lineage_dict_bc' as input. Lineage_dict_bc is a vector of lineage numbers named according to cell barcode. The function outputs the number of true positive predictions (TP), false positive predictions (FP), as well as precision and sensitiivity for various confidence intervals. 

```
> lineage_testing = test_lineage_prediction(lineage_predictions_matrix, lineage_dict_bc)
> lineage_testing
     TP    FP  precision sensitivity
0   274 24062 0.01125904   1.0000000
10  104   128 0.44827586   0.3795620
20   90    82 0.52325581   0.3284672
30   84    56 0.60000000   0.3065693
40   80    36 0.68965517   0.2919708
50   68    22 0.75555556   0.2481752
60   64    14 0.82051282   0.2335766
70   62    10 0.86111111   0.2262774
80   58     2 0.96666667   0.2116788
90   46     0 1.00000000   0.1678832
100  34     0 1.00000000   0.1240876
```

This results can be visualized using `XXX`.

```
visualize_test_result(lineage_testing)
```

<p align="center">
  <img width="500" height="500" src="https://github.com/UPSUTER/GEMLI/Example/GEMIL_GitHub_testing.png">
</p>

## Citation
If you use the package, please cite X
