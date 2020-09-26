options(java.parameters="-Xmx256g")
library(Seurat)
library(scater)
library(readr)
library(dplyr)
library(xlsx)

sample <- c("Org_day1.Seu_mm", "Org_day7.Seu_mm")
dims <- c(6, 8)

for(i in seq(sample)){
  smpl <- strsplit(sample[i], ".Seu_")[[1]][1]
  species <- strsplit(sample[i], ".Seu_")[[1]][2]
 
  # Load gene annotation
  Anno <- read.delim(paste0("data/scRNA/Anno/Anno.", smpl, ".", species, ".biomaRt.txt"), row.names=1)
  
  # Load Seurat Object
  Seu <- read_rds(paste0("data/scRNA/Seurat/RData/", sample[i], ".Epi.SCTransform.Rds"))

  # Run TSNE
  n.dims <- dims[i]
  Seu <- RunTSNE(Seu, dims.use=1:n.dims)

  # tSNE plot for Replicates
  p <- DimPlot(Seu, reduction="tsne", group.by="Dataset")
  pdf(paste0("data/scRNA/Seurat/TSNE/TSNE.", smpl, ".", species, ".Epi.pdf"))
  print(p)
  dev.off()
  
  # Find Clusters
  resolution <- 0.2
  Seu <- FindNeighbors(Seu, dims = 1:n.dims)
  Seu <- FindClusters(Seu, resolution = resolution, verbose=FALSE)
  
  # Save Data
  write_rds(Seu, paste0("data/scRNA/Seurat/RData/", sample[i], ".Epi.SCTransform.SNN.Rds"))

  # Draw tSNE plot
  p <- DimPlot(Seu, reduction="tsne", group.by="seurat_clusters")
  pdf(paste0("data/scRNA/Seurat/TSNE/TSNE.", smpl, ".", species, ".Epi.SCT.SNN.res0.1.pdf"), width=8)
  print(p)
  dev.off()

  # Find Marker Genes (Wilcox)
  markers <- FindAllMarkers(Seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose=FALSE)
  MG <- markers %>% group_by(cluster)
  top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  res <- cbind(Anno[MG$gene, ], as.data.frame(MG))
  res <- res[res$p_val_adj < 0.05,]
  res <- split(res, f=res$cluster)
  
  file_MG <- paste0("data/scRNA/Seurat/MG_list/Marker.Genes.Wilcox.", smpl, ".", species, ".Epi.xlsx")
  for(c in seq(length(res))){
   write.xlsx2(res[[c]], file_MG, row.names=FALSE, append=TRUE, sheetName=paste0("cluster", c-1))
  }

  # Violin Plot for Maker Genes
  dir.create(paste0("data/scRNA/Seurat/VlnPlot/", smpl, ".", species, ".Epi"))
  for(c in seq(length(res))){
    p <- VlnPlot(Seu, features=res[[c]]$gene[1:10], pt.size=0)
    pdf(paste0("data/scRNA/Seurat/VlnPlot/", smpl, ".", species, ".Epi/Cluster", c-1, ".pdf"), width=15, height=12)
    print(p)
    dev.off()
  }
}