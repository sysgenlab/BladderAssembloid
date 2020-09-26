options(java.parameters="-Xmx256g")
library(Seurat)
library(scater)
library(readr)
library(cowplot)
library(dplyr)
library(xlsx)
library(doParallel)

sample <- c("Org_day1.Seu_hs", "Org_day7.Seu_hs")
dims <- c(6, 8)
perplexity <- c(30, 26)

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
  write.table(Anno, file_Anno, sep="\t", quote=F, row.names=FALSE, col.names=TRUE)

  # Load Seurat Object
  Seu <- read_rds(paste0("data/scRNA/Seurat/RData/", sample[i], ".SCTransform.Rds"))

  # Run TSNE
  n.dims <- dims[i]
  Seu <- RunTSNE(Seu, dims.use=1:n.dims, perplexity=perplexity[i])

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
}