options(java.parameters="-Xmx256g")
options(stringsAsFactors=FALSE)
library(Seurat)
library(readr)
library(cowplot)
library(dplyr)
library(xlsx)

Seu <- read_rds("data/scRNA/Seurat/RData/WT-Org.integrated.Seu_mm.SCTransform.Rds")
Anno <- read.delim("data/scRNA/Anno/Anno.Integrated.biomaRt.txt", row.names=1)

# Run TSNE
n.dims <- 9
Seu <- RunTSNE(Seu, dims.use=1:n.dims)

# tSNE plot for Replicates
p <- DimPlot(Seu, reduction="tsne", group.by="Dataset")
pdf(paste0("data/scRNA/Seurat/TSNE/TSNE.WT-Org.mm.pdf"))
print(p)
dev.off()

resolution <- 0.1
Seu <- FindNeighbors(Seu, dims = 1:n.dims)
Seu <- FindClusters(Seu, resolution = resolution)

# Save Data
write_rds(Seu, "data/scRNA/Seurat/RData/WT-Org.Seu_mm.SCTransform.SNN.Rds")

# Draw tSNE plot
p <- DimPlot(Seu, reduction="tsne", group.by="seurat_clusters")
pdf(paste0("data/scRNA/Seurat/TSNE/TSNE.WT-Org.mm.SCT.SNN.res0.1.pdf"), width=8)
print(p)
dev.off()

# Find Marker Genes (Wilcox)
markers <- FindAllMarkers(Seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose=FALSE)
MG <- markers %>% group_by(cluster)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
res <- cbind(Anno[MG$gene, ], as.data.frame(MG))
res <- res[res$p_val_adj < 0.05,]
res <- split(res, f=res$cluster)

file_MG <- paste0("data/scRNA/Seurat/MG_list/Marker.Genes.Wilcox.WT-Org.mm.xlsx")
for(c in seq(length(res))){
  write.xlsx2(res[[c]], file_MG, row.names=FALSE, append=TRUE, sheetName=paste0("cluster", c-1))
}

# Violin Plot for Maker Genes
dir.create(paste0("data/scRNA/Seurat/VlnPlot/WT-Org.mm"))
for(c in seq(length(res))){
  p <- VlnPlot(Seu, features=res[[c]]$gene[1:10], pt.size=0)
  pdf(paste0("data/scRNA/Seurat/VlnPlot/WT-Org.mm/Cluster", c-1, ".pdf"), width=15, height=12)
  print(p)
  dev.off()
}
