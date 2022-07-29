library(Seurat)
library(cowplot)
KO.data <- Read10X(data.dir = "D:/Jmjd1c_Treg_Tumor/single_cell_RNAseq/201692_20210126_1_foxp3jmjd1c_mc205_tumour_wt_ko_treg/Treg_KO/outs/filtered_feature_bc_matrix")
WT.data <- Read10X(data.dir = "D:/Jmjd1c_Treg_Tumor/single_cell_RNAseq/201692_20210126_1_foxp3jmjd1c_mc205_tumour_wt_ko_treg/WT/outs/filtered_feature_bc_matrix")

KO.data <- as.data.frame(KO.data)
WT.data <- as.data.frame(WT.data)

for (i in 1:4105) {
  colnames(KO.data)[i] <- paste(colnames(KO.data)[i],"KO",i,sep = "-")  
}

for (i in 1:3905) {
  colnames(WT.data)[i] <- paste(colnames(WT.data)[i],"WT",i,sep = "-")  
}

KO.metadata<-data.frame(colnames(KO.data),rep("KO",4105))
WT.metadata<-data.frame(colnames(WT.data),rep("WT",3905))
colnames(KO.metadata)<-c("barcode","group")
colnames(WT.metadata)<-c("barcode","group")
rownames(KO.metadata)<-KO.metadata[,1]
rownames(WT.metadata)<-WT.metadata[,1]


KO <- CreateSeuratObject(counts = KO.data, project = "IMMUNE_KO", min.cells = 5,meta.data = KO.metadata)
KO$type <- "KO"
KO <- subset(KO, subset = nFeature_RNA > 500)
KO <- NormalizeData(KO, verbose = FALSE)
KO <- FindVariableFeatures(KO, selection.method = "vst", nfeatures = 2000)

WT <- CreateSeuratObject(counts = WT.data, project = "IMMUNE_WT", min.cells = 5,meta.data = WT.metadata)
WT$type <- "WT"
WT <- subset(WT, subset = nFeature_RNA > 500)
WT <- NormalizeData(WT, verbose = FALSE)
WT <- FindVariableFeatures(WT, selection.method = "vst", nfeatures = 2000)

immune.anchors <- FindIntegrationAnchors(object.list = list(KO, WT), dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)

DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)

immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.8)

p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "type")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
p1
p2
DimPlot(immune.combined, reduction = "umap", split.by = "type")

immune.combined <- RunTSNE(immune.combined, dims = 1:20)
DimPlot(immune.combined, reduction = "tsne",split.by = "type")
DimPlot(immune.combined, reduction = "tsne",group.by  = "type")
DimPlot(immune.combined, reduction = "tsne")
immune.combined.markers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
write.table(immune.combined.markers,file = "D:/Jmjd1c_Treg_Tumor/single_cell_RNAseq/201692_20210126_1_foxp3jmjd1c_mc205_tumour_wt_ko_treg/Seurat_result/immune.combined.markers.txt",sep = "\t")

for (i in 0:3) {
  immune.combined_i.conservemarkers <-FindConservedMarkers(immune.combined, ident.1 = i, grouping.var = "type", verbose = FALSE,min.pct = 0, logfc.threshold = 0)
  write.table(immune.combined_i.conservemarkers,file = paste("D:/Jmjd1c_Treg_Tumor/single_cell_RNAseq/201692_20210126_1_foxp3jmjd1c_mc205_tumour_wt_ko_treg/Seurat_result/immune.combined_conserve",i,"txt",sep = "."),sep = "\t")
}

immune.combined$cell.type <- paste(Idents(immune.combined), immune.combined$type, sep = "_")
immune.combined$cell.cluster <- Idents(immune.combined)
for (i in 0:3) {
  i.cells <- subset(immune.combined, idents = i)
  Idents(i.cells) <- "type"
  avg.i.cells <- log1p(AverageExpression(i.cells, verbose = FALSE)$RNA)
  avg.i.cells$gene <- rownames(avg.i.cells)
  write.table(avg.i.cells,file = paste("D:/Jmjd1c_Treg_Tumor/single_cell_RNAseq/201692_20210126_1_foxp3jmjd1c_mc205_tumour_wt_ko_treg/Seurat_result/avg",i,"txt",sep = "."),sep = "\t")
}
Idents(immune.combined) <- "cell.type"
for (i in 0:3) {
  immune.combined_i.response <- FindMarkers(immune.combined, ident.1 = paste(i,"KO",sep = "_"), ident.2 = paste(i,"WT",sep = "_"), verbose = FALSE,min.pct = 0, logfc.threshold = 0)
  write.table(immune.combined_i.response,file = paste("D:/Jmjd1c_Treg_Tumor/single_cell_RNAseq/201692_20210126_1_foxp3jmjd1c_mc205_tumour_wt_ko_treg/Seurat_result/immune.combined",i,"txt",sep = "."),sep = "\t")
}

immune.combined@meta.data$cell.cluster<-Idents(immune.combined)
Idents(immune.combined)<-immune.combined@meta.data$group
immune.combined.markers<-FindMarkers(immune.combined,ident.1 = "WT",ident.2 = "KO",logfc.threshold = 0,min.pct = 0)
write.table(immune.combined.markers,file = "D:/Jmjd1c_Treg_Tumor/single_cell_RNAseq/201692_20210126_1_foxp3jmjd1c_mc205_tumour_wt_ko_treg/Seurat_result/immune.combined_total.txt",sep = "\t")

library(dplyr)
library(Seurat)
library(patchwork)
top10 <- immune.combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(immune.combined, features = top10$gene) #+ NoLegend()

top100 <- immune.combined.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
DoHeatmap(immune.combined, features = top100$gene) #+ NoLegend()
