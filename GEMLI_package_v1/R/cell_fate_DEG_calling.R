cell_fate_DEG_calling<-function(GEMLI_items, ident1, ident2, min.pct=0.05, logfc.threshold=0.1)
{
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Seurat is required for this function. Please install and load it before using cell_fate_DEG_calling")
    }
    if (class(GEMLI_items)=='list') {
        gene_expression = Matrix::Matrix(GEMLI_items[['gene_expression']])
    } else if (class(GEMLI_items)=='GEMLI') {
        gene_expression = Matrix::Matrix(GEMLI_items@gene_expression)
    } else {
        stop('Object GEMLI_items should be either of class list or GEMLI')
    }
    GEMLI_Seurat<-Seurat::CreateSeuratObject(gene_expression, project = "SeuratProject", assay = "RNA")
    Metadata<-GEMLI_items[['cell_fate_analysis']]
    Metadata$ident<-NA
    Metadata$ident[Metadata$cell.fate %in% ident1]<-"ident1"
    Metadata$ident[Metadata$cell.fate%in%ident2]<-"ident2"
    Meta<-as.data.frame(Metadata[,c(5)])
    rownames(Meta)<-Metadata$cell.ID
    colnames(Meta)<-c("cell.fate")
    GEMLI_Seurat<-Seurat::AddMetaData(GEMLI_Seurat, Meta, col.name = NULL)
    Seurat::DefaultAssay(object = GEMLI_Seurat) <- "RNA"
    Seurat::Idents(GEMLI_Seurat) <- GEMLI_Seurat$cell.fate
    GEMLI_Seurat <- Seurat::NormalizeData(GEMLI_Seurat, normalization.method = "LogNormalize", scale.factor = 10000)
    DEG <- Seurat::FindMarkers(object = GEMLI_Seurat, ident.1 = "ident1", ident.2 = "ident2", min.pct = min.pct, logfc.threshold = logfc.threshold, assay='RNA')
    GEMLI_items[['DEG']]<-DEG
    return(GEMLI_items)
}
