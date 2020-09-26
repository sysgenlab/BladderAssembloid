options(java.parameters="-Xmx256g")
library(Seurat)
library(scater)
library(readr)
library(cowplot)
library(dplyr)
library(xlsx)
library(doParallel)

sample <- c("WT.Seu_mm", "Org_day1.Seu_mm", "Org_day7.Seu_mm")
dims <- c(15, 12, 11)

dir.create("data/scRNA/Anno")
dir.create("data/scRNA/Seurat/TSNE")
dir.create("data/scRNA/Seurat/MG_list")
dir.create("data/scRNA/Seurat/VlnPlot")

for(i in 1:length(sample)){
  smpl <- strsplit(sample[i], ".Seu_")[[1]][1]
  species <- strsplit(sample[i], ".Seu_")[[1]][2]

  # Load & Save Gene Annotation
  sce <- read_rds(paste0("data/scRNA/scater/RData_filt/", smpl, ".sce_", species, ".QC.norm.filt.Rds"))
  rownames(sce) <- rowData(sce)$symbol
  Anno <- rowData(sce)[,1:10]
  Anno$symbol[duplicated(Anno$symbol)] <- paste0(Anno$symbol[duplicated(Anno$symbol)], "_dup")
  rownames(Anno) <- Anno$symbol
  Anno <- Anno[,c(4,5,6,8,9,10)]
  
  file_Anno <- paste0("data/scRNA/Anno/Anno.", smpl, ".", species, ".biomaRt.txt")
  write.table(Anno, file_Anno, sep="\t", quote=F, row.names=TRUE, col.names=TRUE)

  # Load Seurat Object
  Seu <- read_rds(paste0("data/scRNA/Seurat/RData/", sample[i], ".SCTransform.Rds"))

  # Run TSNE
  n.dims <- dims[i]
  Seu <- RunTSNE(Seu, dims.use=1:n.dims)

  # tSNE plot for Replicates
  p <- DimPlot(Seu, reduction="tsne", group.by="Dataset")
  pdf(paste0("data/scRNA/Seurat/TSNE/TSNE.", smpl, ".", species, ".pdf"))
  print(p)
  dev.off()

  # Find Clusters
  resolution <- 0.1
  Seu <- FindNeighbors(Seu, dims = 1:n.dims)
  Seu <- FindClusters(Seu, resolution = resolution)

  # Save Data
  write_rds(Seu, paste0("data/scRNA/Seurat/RData/", sample[i], ".SCTransform.SNN.Rds"))

  # Draw tSNE plot
  p <- DimPlot(Seu, reduction="tsne", group.by="seurat_clusters")
  pdf(paste0("data/scRNA/Seurat/TSNE/TSNE.", smpl, ".", species, ".SCT.SNN.res0.1.pdf"), width=8)
  print(p)
  dev.off()

  # Find Marker Genes (Wilcox)
  if(length(levels(Idents(Seu))) == 1){next}
  markers <- FindAllMarkers(Seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose=FALSE)
  MG <- markers %>% group_by(cluster)
  top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  res <- cbind(Anno[MG$gene, ], MG)
  res <- res[res$p_val_adj < 0.05,]
  res <- split(res, f=res$cluster)
  
  file_MG <- paste0("data/scRNA/Seurat/MG_list/Marker.Genes.Wilcox.", smpl, ".", species, ".xlsx")
  for(c in seq(length(res))){
    write.xlsx2(res[[c]], file_MG, row.names=FALSE, append=TRUE, sheetName=paste0("cluster", c-1))
  }

  # Violin Plot for Maker Genes
  dir.create(paste0("data/scRNA/Seurat/VlnPlot/", smpl, ".", species))
  for(c in seq(length(res))){
    p <- VlnPlot(Seu, features=res[[c]]$gene[1:10], pt.size=0)
    pdf(paste0("data/scRNA/Seurat/VlnPlot/", smpl, ".", species, "/Cluster", c-1, ".pdf"), width=15, height=12)
    print(p)
    dev.off()
  }
}
