# GEMLI: Gene expression memory for lineage identification

GEMLI is an R package to predict cell lineages (cells with a common ancestor) from single cell RNA datasets (or any other type of single-cell gene expression data) and to call genes with a high gene expression memory on the predicted cell lineages. It is described in A.S. Eisele*, M. Tarbier*, A.A. Dormann, V. Pelechano, D.M. Suter | "Barcode-free prediction of cell lineages from scRNA-seq datasets" | 2022 bioRxiv.

The approach is based on findings of Phillips et al. 2019 where it was shown that some genes show varying gene expression across cell lineages that is stable over multiple cell generation.

## Installation
To install the package simply download the folder "GEMLI_package_v0". Open the project "GEMLI.Rproj" in RStudio. Then click in the RStudio menu on "Build" and then "Install and Restart". You can then close the window with the project with `rstudioapi::executeCommand('closeProject')` and use the library in R.

## Development and feedback
We are still working to make GEMLI more intuitive, user-friendly, faster and versatile. Therefore exiting functions will still be updated and new functionalities will be added. We'll publish a list of changes for each version for you to keep track. Your feedback is very welcome and will help us make GEMLI even better. What do you like about GEMLI? Something not working? What functions are you missing? Let us know! Contact us: [marcel.tarbier@scilifelab.se]

## Example 1: small lineages in mouse embryonic stem cells

### Load example data
First we load the example data.

```
> load('GEMLI_example_data_matrix.RData')
> load('GEMLI_example_barcode_information.RData')
```

### Create a GEMLI items list
GEMLI's inputs and outputs are stored in a list of objects with predefined names. To run GEMLI you need at least a quality controlled and normalized gene expression matrix (rows = genes/features, colums = cells/samples). In this example we also provide a ground truth for lineages stemming from a barcoding experiment (values = barcode ID, names = cell IDs).

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
Since we have barcoding data for this dataset we can test the predicted lineages against our ground truth.
The `test_lineage_prediction` function again takes our `GEMLI_items` as input. It's important that a predcition has been run first with `predict_lineages`. It outputs the number of true positive predictions (TP), false positive predictions (FP), as well as precision and sensitivity for various confidence intervals. The output can be visualized by setting `plot_results` to `true`/`T`.

```
> GEMLI_items = test_lineages(GEMLI_items)
>
> GEMLI_items$testing_results
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
> test_lineages(GEMLI_items, plot_results=T)
```

<p align="center">
  <img width="500" height="500" src="https://github.com/UPSUTER/GEMLI/blob/main/Example/GEMIL_GitHub_testing.png">
</p>

### Visualize predictions as network
We can also investigate our predictions by visualizing them as a network with the `visualize_as_network` function. Here we need to set a `cutoff` that defines which predictions we want to consider. It represents a confidence score and high values yield fewer predictions with high precision while low values yield more predcitions with lower precision.

```
> visualize_as_network(GEMLI_items, cutoff=90) # left image
> visualize_as_network(GEMLI_items, cutoff=50) # right image
```
<p float="left">
  <img width="500" height="500" src="https://github.com/UPSUTER/GEMLI/blob/main/Example/GEMLI_GitHub_network_90.png">
  <img width="500" height="500" src="https://github.com/UPSUTER/GEMLI/blob/main/Example/GEMLI_GitHub_network_50.png">
</p>

If a ground truth e.g. from barcoding is avalable we can set `ground_truth` to `true`/`T` to highlight false predictions with red edges. Cells without barcode information will be displaye in white.

```
> visualize_as_network(GEMLI_items, cutoff=90, ground_truth=T) # left image
> visualize_as_network(GEMLI_items, cutoff=50, ground_truth=T) # right image
```
<p float="left">
  <img width="500" height="500" src="https://github.com/UPSUTER/GEMLI/blob/main/Example/GEMLI_GitHub_network_90_GT.png">
  <img width="500" height="500" src="https://github.com/UPSUTER/GEMLI/blob/main/Example/GEMLI_GitHub_network_50_GT.png">
</p>

### Extract lineage information
Now we can extract the lineage information with the `prediction_to_lineage_information` function. Again we need to set a `cutoff` that defines which predictions we want to consider. The function outputs both a lineage table and a 'dictionary', a vector that has the lineage number as values and the cell IDs as names.

```
> GEMLI_items = prediction_to_lineage_information(GEMLI_items, cutoff=50)
>
> GEMLI_items$predicted_lineage_table[1:5,]
     cell.ID              clone.ID
[1,] "AAACGAACAGGTGTGA-1" "1"
[2,] "AAAGGTAGTTGCTTGA-1" "2"
[3,] "AACCACAAGTTTGTCG-1" "3"
[4,] "AAGCCATGTTCCACGG-1" "4"
[5,] "AAGCGAGGTACGGCAA-1" "5"
>
> GEMLI_items$predicted_lineages[1:5,]
AAACGAACAGGTGTGA-1 AAAGGTAGTTGCTTGA-1 AACCACAAGTTTGTCG-1 AAGCCATGTTCCACGG-1 AAGCGAGGTACGGCAA-1
                 1                  2                  3                  4                  5
```

### Trim lineages that are too big

In some applications it may be useful to trim lineages that are too big. For instance if it is known that cells should have undergone only a certain number of divisions effectively limiting the lineage size or if you are only interested in sister cell pairs. Similarly, if you investigate large lineages you want to avoid lineages being merged due to few false predictions between otherwise well-interconnected lineages. The `suggest_network_trimming_to_size` function allows you to preview what a trimming to size would look like. It will again show the predicted lineages as networks but highlight all connections that would be trimmed given a certain size restriction (`max_size`). If you are happy with the suggested trimming you create new trimmed `GEMLI_items` list, in this example we called it `GEMLI_items_post_processed`. You can then again visualize the predictions to see how the changes have affected your predictions.
```
> suggest_network_trimming_to_size(GEMLI_items, max_size=2, cutoff=50) # left image
> GEMLI_items_post_processed = trim_network_to_size(GEMLI_items, max_size=2, cutoff=50)
> visualize_as_network(GEMLI_items_post_processed, cutoff=50) # right image
```
<p float="left">
  <img width="500" height="500" src="https://github.com/UPSUTER/GEMLI/blob/main/Example/GEMLI_GitHub_network_50_GT_ST.png">
  <img width="500" height="500" src="https://github.com/UPSUTER/GEMLI/blob/main/Example/GEMLI_GitHub_network_50_trim.png">
</p>

## Citation
If you use the package, please cite A.S. Eisele*, M. Tarbier*, A.A. Dormann, V. Pelechano, D.M. Suter | "Barcode-free prediction of cell lineages from scRNA-seq datasets" | 2022 bioRxiv.
