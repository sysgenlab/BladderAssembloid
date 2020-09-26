options(java.parameters="-Xmx32g")
library(Seurat)
library(scater)
library(readr)
library(dplyr)
library(xlsx)
library(cowplot)

sample <- c("Org_day1.Seu_mm", "Org_day7.Seu_mm")

for(i in seq(sample)){
  smpl <- strsplit(sample[i], ".Seu_")[[1]][1]
  species <- strsplit(sample[i], ".Seu_")[[1]][2]

  sce <- read_rds(paste0("data/scRNA/scater/RData_filt/", smpl, ".sce_", species, ".QC.norm.filt.Rds"))

  Seu <- as.Seurat(sce)
  colnames(Seu@meta.data)[which(colnames(Seu@meta.data) == "total_counts")] <- "nCount_RNA"
  colnames(Seu@meta.data)[which(colnames(Seu@meta.data) == "total_features_by_counts")] <- "nFeature_RNA"
  Seu_type <- read_rds(paste0("data/scRNA/Seurat/RData/", sample[i], ".SCTransform.SNN.Rds"))

  # Epithelial cells population
  Seu <- Seu[,Seu_type$Type == "Epithelial cells"]

  # Calculate Mitochonrial gene proportions
  Seu <- PercentageFeatureSet(Seu, pattern = "^mt-", col.name="percent.mt")

  # SCTransform
  Seu <- SCTransform(Seu, vars.to.regress = "percent.mt")

  # PCA
  Seu <- RunPCA(object = Seu, verbose=FALSE)

  # Exploring Optimal Clustering Number
  Seu <- JackStraw(Seu, num.replicate=100, verbose=FALSE)
  Seu <- ScoreJackStraw(Seu, dims = 1:20)
  
  pdf(paste0("data/scRNA/Seurat/Dimensionality/", smpl, ".", species, ".Epi.JackStrawPlot.ElbowPlot.pdf"), height=7, width=16)
    print(plot_grid(JackStrawPlot(Seu, dims = 1:20), ElbowPlot(Seu)))
  dev.off()

  write_rds(Seu, paste0("data/scRNA/Seurat/RData/", sample[i], ".Epi.SCTransform.Rds"))
}