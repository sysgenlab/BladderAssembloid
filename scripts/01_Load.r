# The script follows previous study: http://oshlacklab.com/combes-organoid-paper/

### Quality Control ###
# Read in single-cell RNA-Seq data, produce various quality control plots and remove any low-quality cells or uniformative genes.
## The number of counts per barcode
## The number of genes per barcode
## The fraction of counts from mitochondrial genes per barcode

library(scater);library(DropletUtils);library(scran)  # scRNA-seq
library(stringr)
library(biomaRt)
library(cowplot)
library(BiocParallel)  # Parallel
library(BiocGenerics)
library(jsonlite)  # Output
library(tidyverse)  # Tidyverse

dir.create("data/scRNA/scater")
dir.create("data/scRNA/scater/RData")
dir.create("data/scRNA/scater/tSNE")

sample <- c("WT", "Org_day1", "Org_day7")

for(i in 1:length(sample)){
  # Read in Single-Cell RNA-Seq data
  ##DropletUtils-1.4.3 - read10xCounts
  ##stringr-1.4.0 - str_sub
  ##SummarizedExperiment-1.14.1 - colData, rowData
  sce <- read10xCounts(paste0("data/scRNA/matrix/", sample[i], "/outs/filtered_feature_bc_matrix"), col.names=TRUE)
  class <- read.csv(paste0("data/scRNA/matrix/", sample[i], "/outs/analysis/gem_classification.csv"))
  cell.names <- paste0("Cell", 1:ncol(sce))
  colnames(sce) <- cell.names
  rownames(class) <- cell.names
  barcodes <- colData(sce)$Barcode
  colData(sce)$Cell <- cell.names
  colData(sce)$Barcode <- str_sub(barcodes, end = -3)
  colData(sce)$Sample <- str_sub(barcodes, start = -1)
  colData(sce)$Dataset <- sample[i]
  colnames(rowData(sce)) <- c("gene_id", "symbol", "type")

  multiplet <- (class$GRCh38+1)/(class$mm10+1) < 3/2 & (class$GRCh38+1)/(class$mm10+1) > 2/3
  class <- class[!multiplet,]
  sce <- sce[,rownames(class)]
  sce_hs <- sce[grep("^GRCh38", rownames(sce)), rownames(class[class$GRCh38 > class$mm10,])]
  sce_mm <- sce[grep("^mm10", rownames(sce)), rownames(class[class$mm10 > class$GRCh38,])]
  class_hs <- class[colnames(sce_hs), ]
  class_mm <- class[colnames(sce_mm), ]

  rownames(sce_hs) <- sub("^GRCh38_", "", rownames(sce_hs))
  rowData(sce_hs)$gene_id <- sub("^GRCh38_", "", rowData(sce_hs)$gene_id)
  rowData(sce_hs)$symbol <- sub("^GRCh38_", "", rowData(sce_hs)$symbol)
  colData(sce_hs)$call <- class_hs$call
  rownames(sce_mm) <- sub("^mm10___", "", rownames(sce_mm))
  rowData(sce_mm)$gene_id <- sub("^mm10___", "", rowData(sce_mm)$gene_id)
  rowData(sce_mm)$symbol <- sub("^mm10___", "", rowData(sce_mm)$symbol)
  colData(sce_mm)$call <- class_mm$call

  # Annotation
  ## Adding Features from biomaRt
  ##biomaRt-2.40.4 - useMart, getBM
#  org <- "mmusculus_gene_ensembl"
  id <- "ensembl_gene_id"
#  attr <- c("ensembl_gene_id", "entrezgene_id", "external_gene_name", "mgi_symbol", "chromosome_name", "description", "gene_biotype", "percentage_gene_gc_content")
  values <- rownames(sce_mm)
#  bmart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset=org, host="www.ensembl.org")
#  coords_mm <- getBM(attributes=attr, filters=id, values = values, mart=bmart)
#  write_rds(coords_mm, "data/scRNA/scater/RData/coords.biomaRt.mm.Rds")
  coords_mm <- read_rds("data/scRNA/scater/RData/coords.biomaRt.mm.Rds")
  
  mm <- match(values, coords_mm[[id]])
  coords_all <- coords_mm[mm, ]
  old_rdata <- rowData(sce_mm)
  keep <- !(colnames(old_rdata) %in% colnames(coords_all))
  new_rdata <- cbind(old_rdata[, keep], coords_all)
  new_rdata$description <- gsub("\\s\\[.*\\]", "", new_rdata$description)
  rowData(sce_mm) <- new_rdata
  
#  org <- "hsapiens_gene_ensembl"
#  attr <- c("ensembl_gene_id", "entrezgene_id", "external_gene_name", "hgnc_symbol", "chromosome_name", "description", "gene_biotype", "percentage_gene_gc_content")
  values <- rownames(sce_hs)
#  bmart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset=org, host="www.ensembl.org")
#  coords_hs <- getBM(attributes=attr, filters=id, values = values, mart=bmart)
#  write_rds(coords_hs, "data/scRNA/scater/RData/coords.biomaRt.hs.Rds")
  coords_hs <- read_rds("data/scRNA/scater/RData/coords.biomaRt.hs.Rds")

  mm <- match(values, coords_hs[[id]])
  coords_all <- coords_hs[mm, ]
  old_rdata <- rowData(sce_hs)
  keep <- !(colnames(old_rdata) %in% colnames(coords_all))
  new_rdata <- cbind(old_rdata[, keep], coords_all)
  new_rdata$description <- gsub("\\s\\[.*\\]", "", new_rdata$description)
  rowData(sce_hs) <- new_rdata


  ## Adding Cell Cycle Phase
  ##scran-1.12.1 - cyclone
  cc.pairs_mm <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package = "scran"))
  cc.pairs_hs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package = "scran"))
  cycles_mm <- cyclone(sce_mm, pairs = cc.pairs_mm, BPPARAM = MulticoreParam(workers=30), verbose = TRUE)
  cycles_hs <- cyclone(sce_hs, pairs = cc.pairs_hs, BPPARAM = MulticoreParam(workers=30), verbose = TRUE)

  colData(sce_mm)$G1Score <- cycles_mm$scores$G1
  colData(sce_mm)$SScore <- cycles_mm$scores$G1
  colData(sce_mm)$G2MScore <- cycles_mm$scores$G2M
  colData(sce_mm)$CellCycle <- cycles_mm$phases

  colData(sce_hs)$G1Score <- cycles_hs$scores$G1
  colData(sce_hs)$SScore <- cycles_hs$scores$G1
  colData(sce_hs)$G2MScore <- cycles_hs$scores$G2M
  colData(sce_hs)$CellCycle <- cycles_hs$phases  


  ## Calculating CPM
  ##scater-1.12.2 - calculateCPM
  cpm(sce_mm) <- calculateCPM(sce_mm, use_size_factors = FALSE)
  cpm(sce_hs) <- calculateCPM(sce_hs, use_size_factors = FALSE)


  ## Calculating QC metrics
  ##scater-1.12.2 - calculateQCMetrics
  sce_mm <- calculateQCMetrics(sce_mm, BPPARAM=MulticoreParam(workers=30))
  colData(sce_mm)$pct_dropout <- 100 * (1 - colData(sce_mm)$total_features_by_counts / nrow(sce_mm))
  sce_hs <- calculateQCMetrics(sce_hs, BPPARAM=MulticoreParam(workers=30))
  colData(sce_hs)$pct_dropout <- 100 * (1 - colData(sce_hs)$total_features_by_counts / nrow(sce_hs))


  ## Adding Percentage of Mitochondrial Genes
  ##BiocGenerics-0.30.0 - counts
  mt <- rowData(sce_mm)$chromosome_name == "MT"
  mt[is.na(mt)] <- FALSE
  mt <- mt | grepl("mitochondrial", rowData(sce_mm)$description)
  colData(sce_mm)$PctCountsMT <- colSums(as.matrix(counts(sce_mm)[mt, ])) / colData(sce_mm)$total_counts * 100

  mt <- rowData(sce_hs)$chromosome_name == "MT"
  mt[is.na(mt)] <- FALSE
  mt <- mt | grepl("mitochondrial", rowData(sce_hs)$description)
  colData(sce_hs)$PctCountsMT <- colSums(as.matrix(counts(sce_hs)[mt, ])) / colData(sce_hs)$total_counts * 100


  ## Adding Percentage of Ribosomal RNAs
  is.ribo <- grepl("ribosom", rowData(sce_mm)$description)
  colData(sce_mm)$PctCountsRibo <- colSums(as.matrix(counts(sce_mm)[is.ribo, ])) / colData(sce_mm)$total_counts * 100

  is.ribo <- grepl("ribosom", rowData(sce_hs)$description)
  colData(sce_hs)$PctCountsRibo <- colSums(as.matrix(counts(sce_hs)[is.ribo, ])) / colData(sce_hs)$total_counts * 100

  ## Saving Results
  write_rds(sce_mm, paste0("data/scRNA/scater/RData/", sample[i], ".sce_mm.QC.Rds"))
  write_rds(sce_hs, paste0("data/scRNA/scater/RData/", sample[i], ".sce_hs.QC.Rds"))

  ## Normalization
  sce_mm <- scater::normalize(sce_mm)
  sce_mm <- runPCA(sce_mm)
  sce_mm <- runTSNE(sce_mm)

  sce_hs <- scater::normalize(sce_hs)
  sce_hs <- runPCA(sce_hs)
  sce_hs <- runTSNE(sce_hs)

  write_rds(sce_mm, paste0("data/scRNA//RData/", sample[i], ".sce_mm.QC.norm.Rds"))
  write_rds(sce_hs, paste0("data/scRNA/scater/RData/", sample[i], ".sce_hs.QC.norm.Rds"))
  
  p1 <- plotTSNE(sce_mm, colour_by = "log10_total_counts") +
    ggtitle("Total Counts")
  p2 <- plotTSNE(sce_mm, colour_by = "CellCycle") +
    ggtitle("Cell Cycle")
  p3 <- plotTSNE(sce_mm, colour_by = "pct_dropout") +
    ggtitle("Dropout")
  p4 <- plotTSNE(sce_mm, colour_by = "PctCountsMT") +
    ggtitle("Mitochondrial Genes")
  p5 <- plotTSNE(sce_mm, colour_by = "PctCountsRibo") +
    ggtitle("Ribosomal Genes")
  
  pdf(paste0("data/scRNA/scater/tSNE/", sample[i], ".QC.mm.tSNE.pdf"), width=10, height=20)
  print(plot_grid(p1, p2, p3, p4, p5, ncol = 2))
  dev.off()

  p1 <- plotTSNE(sce_hs, colour_by = "log10_total_counts") +
    ggtitle("Total Counts")
  p2 <- plotTSNE(sce_hs, colour_by = "CellCycle") +
    ggtitle("Cell Cycle")
  p3 <- plotTSNE(sce_hs, colour_by = "pct_dropout") +
    ggtitle("Dropout")
  p4 <- plotTSNE(sce_hs, colour_by = "PctCountsMT") +
    ggtitle("Mitochondrial Genes")
  p5 <- plotTSNE(sce_hs, colour_by = "PctCountsRibo") +
    ggtitle("Ribosomal Genes")

  pdf(paste0("data/scRNA/scater/tSNE/", sample[i], ".QC.hs.tSNE.pdf"), width=10, height=20)
  print(plot_grid(p1, p2, p3, p4, p5, ncol = 2))
  dev.off()
}
