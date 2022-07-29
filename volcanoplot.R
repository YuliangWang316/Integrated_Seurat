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

immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:4)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:4)
immune.combined <- FindClusters(immune.combined, resolution = 0.1)

p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "type")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
p1
p2
DimPlot(immune.combined, reduction = "umap", split.by = "type")

immune.combined <- RunTSNE(immune.combined, dims = 1:4)
DimPlot(immune.combined, reduction = "tsne",split.by = "type",label = TRUE)
DimPlot(immune.combined, reduction = "tsne",group.by  = "type")
DimPlot(immune.combined, reduction = "tsne", label = TRUE)

immune.combined@meta.data$cell.cluster<-Idents(immune.combined)
Idents(immune.combined)<-immune.combined@meta.data$group
immune.combined.markers<-FindMarkers(immune.combined,ident.1 = "KO",ident.2 = "WT",logfc.threshold = 0,min.pct = 0)
write.table(immune.combined.markers,file = "D:/Jmjd1c_Treg_Tumor/single_cell_RNAseq/201692_20210126_1_foxp3jmjd1c_mc205_tumour_wt_ko_treg/Seurat_result/immune.combined_total.txt",sep = "\t")

library(dplyr)
library(Seurat)
library(patchwork)
top10 <- immune.combined.markers  %>% top_n(n = 10, wt = avg_logFC)
top10$gene<-rownames(top10)
DoHeatmap(immune.combined, features = top10$gene) #+ NoLegend()

top2000 <- immune.combined.markers  %>% top_n(n = 2000, wt = avg_logFC)
top2000$gene<-rownames(top2000)
DoHeatmap(immune.combined, features = top2000$gene) #+ NoLegend()
a<-immune.combined@assays[["integrated"]]@scale.data
write.table(a,file = "d:/Jmjd1c_Treg_Tumor/single_cell_RNAseq/201692_20210126_1_foxp3jmjd1c_mc205_tumour_wt_ko_treg/Seurat_result/intescaledata.txt",sep = "\t")

gene_list=immune.combined.markers[,c("avg_logFC","p_val_adj")]
#gene_list=log10(gene_list)
#gene_list[,"log2FoldChange"]=log2(gene_list[,"log2FoldChange"])
colnames(gene_list)=c("logFC","padj")
gene_list$change = ifelse(gene_list$padj < 0.000001 & abs(gene_list$logFC) >= 0.1, 
                             ifelse(gene_list$logFC> 0.1 ,'Up 146 genes','Down 15 genes'),
                             'no_change')
#colored_point<-gene_list[gene_list$threshold == "TRUE",]
library("ggplot2")
pdf("d:/Jmjd1c_Treg_Tumor/single_cell_RNAseq/201692_20210126_1_foxp3jmjd1c_mc205_tumour_wt_ko_treg/Seurat_result/volcanoplot.pdf")
g = ggplot(data=gene_list, aes(x=logFC, y=-log10(padj),color=change)) + geom_point(alpha=0.4, size=1.75)  + xlim(c(-0.5, 1)) + ylim(c(0, 310)) +xlab("log2 fold change ko_vs_wt") + ylab("-log10 p-value") + scale_color_manual(values=c("#191970","#000000","#B0171F")) + theme_set(theme_bw()) + theme(panel.grid.major=element_line(colour=NA))
#+geom_text(mapping=aes(label=rownames(colored_point)),data = colored_point,check_overlap = TRUE)
print(g)
dev.off()
