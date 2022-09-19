# GEMLI: Gene expression memory for lineage identification

GEMLI is an R package to predict cell lineages (cells with a common ancestor) from single cell RNA datasets (or any other type of single-cell gene expression data) and to call genes with a high gene expression memory on the predicted cell lineages. It is described in A.S. Eisele*, M. Tarbier*, A.A. Dormann, V. Pelechano, D.M. Suter | "Barcode-free prediction of cell lineages from scRNA-seq datasets" | 2022 bioRxiv.

The approach is based on findings of Phillips et al. 2019 where it was shown that some genes show varying gene expression across cell lineages that is stable over multiple cell generation.

## Installation
To install the package simply download the folder "Package_test". Open the project "LineAGED.Rproj" in RStudio. Then click in the RStudio menu on "Build" and then "Install and Restart". You can then close the window with the project and use the library in R.

## Getting started

### Load example data
First we load the example data.

```
> load('GEMLI_example_data_matrix.RData')
> load('GEMLI_example_barcode_information.RData')
```

### Create a GEMLI items list
GEMLI's inputs and outputs are stored in a list of objects with predifined names. To run GEMLI you need at least a quality controlled and normalized gene expression matrix (rows = genes/features, colums = cells/samples). In this example we also provide a ground truth for lineages stemming from a barcoding experiment (values = barcode ID, colums = cell IDs).

```
> GEMLI_items = list()
> GEMLI_items[['gene_expression']] = data_matrix
> GEMLI_items[['barcodes']] = lineage_dict_bc
>
> GEMLI_items[['gene_expression']][9:14,1:5]
                   AAACGAACAGGTGTGA-1 AAAGGTAGTTGCTTGA-1 AACCACAAGTTTGTCG-1 AAGCCATGTTCCACGG-1 AAGCGAGGTACGGCAA-1
ENSMUSG00000033845          14.761746         12.9570026          13.240645          8.8596794          12.791617
ENSMUSG00000025903           3.163231          1.2340002           4.878132          0.7383066           0.000000
ENSMUSG00000033813           6.326463          5.5530011           4.181256          8.8596794           5.482122
ENSMUSG00000002459           0.000000          0.0000000           0.000000          0.0000000           0.000000
ENSMUSG00000085623           0.000000          0.0000000           0.000000          0.0000000           0.000000
ENSMUSG00000033793           3.163231          0.6170001           1.393752          2.2149199           5.482122
>
> GEMLI_items[['barcodes']][1:5]
CACAGATAGTGATGGC-1 TATCTTGGTACGGGAT-1 AAACGAACAGGTGTGA-1 AGAGAATAGGTCATAA-1 GAGTGAGTCCAGTACA-1
                 2                  2                  2                  7                  7
```

### Perform lineage prediction
We can then identify cell lineages through repeated iterative clustering (this may take 2-3min). The `predict_lineages` function takes our GEMLI_items as input. It outputs a matrix of all cells against all cells with values corresponding to a confidence score that they are part of the same lineage. 

```
> lineage_predictions_matrix = predict_lineages(GEMLI_items)
>
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
> lineage_testing = test_lineage_prediction(GEMLI_items)
>
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
>
> lineage_testing = test_lineage_prediction(lineage_predictions_matrix, lineage_dict_bc)
```

<p align="center">
  <img width="500" height="500" src="https://github.com/UPSUTER/GEMLI/blob/main/Example/GEMIL_GitHub_testing.png">
</p>

## Citation
If you use the package, please cite A.S. Eisele*, M. Tarbier*, A.A. Dormann, V. Pelechano, D.M. Suter | "Barcode-free prediction of cell lineages from scRNA-seq datasets" | 2022 bioRxiv.
