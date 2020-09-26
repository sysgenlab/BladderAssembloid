# tSNE after SCTransform
## Seurat-3.1.0
options(java.parameters="-Xmx32g")
library(Seurat)
library(scater)
library(readr)
library(cowplot)
library(dplyr)

dir.create("data/scRNA/Seurat")
dir.create("data/scRNA/Seurat/RData")
dir.create("data/scRNA/Seurat/Dimensionality")

sample <- c("WT", "Org_day1", "Org_day7")

## For Mouse data
for(i in 1:length(sample)){
  # Load Data
  sce <- read_rds(paste0("data/scRNA/scater/RData_filt/", sample[i], ".sce_mm.QC.norm.filt.Rds"))
  rownames(sce) <- rowData(sce)$symbol
  Seu <- as.Seurat(sce)
  Anno <- rowData(sce)[,1:10]
  Anno$symbol[duplicated(Anno$symbol)] <- paste0(Anno$symbol[duplicated(Anno$symbol)], "_dup")
  rownames(Anno) <- Anno$symbol

  colnames(Seu@meta.data)[which(colnames(Seu@meta.data) == "total_counts")] <- "nCount_RNA"
  colnames(Seu@meta.data)[which(colnames(Seu@meta.data) == "total_features_by_counts")] <- "nFeature_RNA"

  # Calculate Mitochonrial gene proportions
  Seu <- PercentageFeatureSet(Seu, pattern = "^mt-", col.name="percent.mt")

  # SCTransform
  Seu <- SCTransform(Seu, vars.to.regress = "percent.mt")

  # PCA
  Seu <- RunPCA(object = Seu, verbose=FALSE)

  # Exploring Optimal Clustering Number
  Seu <- JackStraw(Seu, num.replicate=100, verbose=FALSE)
  Seu <- ScoreJackStraw(Seu, dims = 1:20)

  p <- print(plot_grid(JackStrawPlot(Seu, dims = 1:20), ElbowPlot(Seu)))
  pdf(paste0("data/scRNA/Seurat/Dimensionality/", sample[i], ".mm.JackStrawPlot.ElbowPlot.pdf"), height=7, width=16)
  print(p)
  dev.off()

  write_rds(Seu, paste0("data/scRNA/Seurat/RData/", sample[i], ".Seu_mm.SCTransform.Rds"))
}

## For Human data
for(i in grep("Org", sample)){
  # Load Data
  sce <- read_rds(paste0("data/scRNA/scater/RData_filt/", sample[i], ".sce_hs.QC.norm.filt.Rds"))
  rownames(sce) <- rowData(sce)$symbol
  Seu <- as.Seurat(sce)
  Anno <- rowData(sce)[,1:10]
  Anno$symbol[duplicated(Anno$symbol)] <- paste0(Anno$symbol[duplicated(Anno$symbol)], "_dup")
  rownames(Anno) <- Anno$symbol

  colnames(Seu@meta.data)[which(colnames(Seu@meta.data) == "total_counts")] <- "nCount_RNA"
  colnames(Seu@meta.data)[which(colnames(Seu@meta.data) == "total_features_by_counts")] <- "nFeature_RNA"

  # Calculate Mitochonrial gene proportions
  Seu <- PercentageFeatureSet(Seu, pattern = "^MT-", col.name="percent.mt")

  # SCTransform
  Seu <- SCTransform(Seu, vars.to.regress = "percent.mt")

  # PCA
  Seu <- RunPCA(object = Seu, verbose=FALSE)

  # Exploring Optimal Clustering Number
  Seu <- JackStraw(Seu, num.replicate=100, verbose=FALSE)
  Seu <- ScoreJackStraw(Seu, dims = 1:20)

  pdf(paste0("data/scRNA/Seurat/Dimensionality/", sample[i], ".hs.JackStrawPlot.ElbowPlot.pdf"), height=7, width=16)
    print(plot_grid(JackStrawPlot(Seu, dims = 1:20), ElbowPlot(Seu)))
  dev.off()

  write_rds(Seu, paste0("data/scRNA/Seurat/RData/", sample[i], ".Seu_hs.SCTransform.Rds"))
}
