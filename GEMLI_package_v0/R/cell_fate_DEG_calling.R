cell_fate_DEG_calling<-function(GEMLI_items, ident1, ident2, min.pct=0.05, logfc.threshold=0.1)
{
  GEMLI_Seurat<-CreateSeuratObject(GEMLI_items[['gene_expression']], project = "SeuratProject", assay = "RNA")
  Metadata<-GEMLI_items[['cell_fate_analysis']];Metadata$ident<-NA
  Metadata$ident[Metadata$cell.fate %in% ident1]<-"ident1"
  Metadata$ident[Metadata$cell.fate%in%ident2]<-"ident2"
  Meta<-as.data.frame(Metadata[,c(5)]); rownames(Meta)<-Metadata$cell.ID; colnames(Meta)<-c("cell.fate")
  GEMLI_Seurat<-AddMetaData(GEMLI_Seurat, Meta, col.name = NULL); DefaultAssay(object = GEMLI_Seurat) <- "RNA"; Idents(GEMLI_Seurat) <- GEMLI_Seurat$cell.fate
  DEG <- FindMarkers(object = GEMLI_Seurat, ident.1 = "ident1", ident.2 = "ident2", min.pct =min.pct, logfc.threshold = logfc.threshold)
  GEMLI_items[['DEG']]<-DEG
  return(GEMLI_items)
}
