# The script follows previous study: http://oshlacklab.com/combes-organoid-paper/

### Quality Control ###
# Read in single-cell RNA-Seq data, produce various quality control plots and remove any low-quality cells or uniformative genes.
## The number of counts per barcode
## The number of genes per barcode
## The fraction of counts from mitochondrial genes per barcode

library(readr)
library(scater)
library(cowplot)

dir.create("data/scRNA/scater/QCplot")

sample <- c("WT", "Org_day1", "Org_day7")

for(i in 1:length(sample)){
  sce_mm <- read_rds(paste0("data/scRNA/scater/RData/", sample[i], ".sce_mm.QC.norm.Rds"))
  sce_hs <- read_rds(paste0("data/scRNA/scater/RData/", sample[i], ".sce_hs.QC.norm.Rds"))

  colData(sce_mm)$Filtered <- FALSE
  colData(sce_hs)$Filtered <- FALSE

  dir.create(paste0("data/scRNA/scater/QCplot/", sample[i]))

  # Total Counts: Depth
  ##scater-1.12.2 - plotColData
  pdf(paste0("data/scRNA/scater/QCplot/", sample[i], "/01_Total_Counts.mm.pdf"), height=7, width=5)
  p <- plotColData(sce_mm, x = "Sample", y = "total_counts") +
        labs(x="Mouse", y="Count depth") + scale_x_discrete(labels = sample[i]) +
        theme(axis.text=element_text(size=13), axis.title=element_text(size=15))
  print(p)
  dev.off()

  pdf(paste0("data/scRNA/scater/QCplot/", sample[i], "/01_Total_Counts.hs.pdf"), height=7, width=5)
  p <- plotColData(sce_hs, x = "Sample", y = "total_counts") +
        labs(x="Human", y="Count depth") + scale_x_discrete(labels = sample[i]) +
        theme(axis.text=element_text(size=13), axis.title=element_text(size=15)) 
  print(p)
  dev.off()
  
  # Number of Features by Library Size
  ##scater-1.12.2 - plotColData
  ## Cells that express many genes are potential multiplets
  pdf(paste0("data/scRNA/scater/QCplot/", sample[i], "/02_Features_by_Library_Size.mm.pdf"), height=7, width=9)
  p <- plotColData(sce_mm, x = "total_counts", y = "total_features_by_counts", colour_by = "pct_dropout") + 
        labs(x="Count depth", y="Number of genes") +
        theme(axis.text=element_text(size=13), axis.title=element_text(size=15))
  print(p)
  dev.off()

  pdf(paste0("data/scRNA/scater/QCplot/", sample[i], "/02_Features_by_Library_Size.hs.pdf"), height=7, width=9)
  p <- plotColData(sce_hs, x = "total_counts", y = "total_features_by_counts", colour_by = "pct_dropout") + 
        labs(x="Count depth", y="Number of genes") +
        theme(axis.text=element_text(size=13), axis.title=element_text(size=15))
  print(p)
  dev.off()
  
  # Mitochondrial genes
  ## Overexpression of mitochondrial genes can be an indication that cell is stressed or damaged in some way.
  pdf(paste0("data/scRNA/scater/QCplot/", sample[i], "/03_Mitochondrial_Genes.mm.pdf"), height=7, width=5)
  p <- plotColData(sce_mm, x = "Sample", y = "PctCountsMT") + 
        labs(x="Mouse", y="% counts MT") + scale_x_discrete(labels = sample[i]) +
        theme(axis.text=element_text(size=13), axis.title=element_text(size=15))
  print(p)
  dev.off()

  pdf(paste0("data/scRNA/scater/QCplot/", sample[i], "/03_Mitochondrial_Genes.hs.pdf"), height=7, width=5)
  p <- plotColData(sce_hs, x = "Sample", y = "PctCountsMT") +               
        labs(x="Human", y="% counts MT") + scale_x_discrete(labels = sample[i]) +
        theme(axis.text=element_text(size=13), axis.title=element_text(size=15))
  print(p)
  dev.off()
  
  
  # Ribosomal genes
  pdf(paste0("data/scRNA/scater/QCplot/", sample[i], "/04_Ribosomal_Genes.mm.pdf"), height=7, width=5)
  p <- plotColData(sce_mm, x = "Sample", y = "PctCountsRibo") + 
        labs(x="Mouse", y="% counts risomal") + scale_x_discrete(labels = sample[i]) +
        theme(axis.text=element_text(size=13), axis.title=element_text(size=15))
  print(p)
  dev.off()

  pdf(paste0("data/scRNA/scater/QCplot/", sample[i], "/04_Ribosomal_Genes.hs.pdf"), height=7, width=5)
  p <- plotColData(sce_hs, x = "Sample", y = "PctCountsRibo") +         
        labs(x="Human", y="% counts risomal") + scale_x_discrete(labels = sample[i]) +
        theme(axis.text=element_text(size=13), axis.title=element_text(size=15))
  print(p)
  dev.off()
 
    
  # Housekeeping genes
  actb.id <- rowData(sce_mm)[which(rowData(sce_mm)$symbol == "Actb"),1]
  gapdh.id <- rowData(sce_mm)[which(rowData(sce_mm)$symbol == "Gapdh"),1]
  key <- c("Actb", "Gapdh")
  names(key) <- c(actb.id, gapdh.id)
  
  pdf(paste0("data/scRNA/scater/QCplot/", sample[i], "/05_Housekeeping_Genes.mm.pdf"), height=7, width=8)
  p <- plotExpression(sce_mm, gapdh.id, x = actb.id) +
        labs(x="Actb", y="Gapdh") +
        theme(axis.text=element_text(size=13), axis.title=element_text(size=15),
            strip.background = element_blank(), strip.text.x = element_blank())
  print(p)
  dev.off()

  actb.id <- rowData(sce_hs)[which(rowData(sce_hs)$symbol == "ACTB"),1]
  gapdh.id <- rowData(sce_hs)[which(rowData(sce_hs)$symbol == "GAPDH"),1]
  key <- c("ACTB", "GAPDH")
  names(key) <- c(actb.id, gapdh.id)

  pdf(paste0("data/scRNA/scater/QCplot/", sample[i], "/05_Housekeeping_Genes.hs.pdf"), height=7, width=8)
  p <- plotExpression(sce_hs, gapdh.id, x = actb.id) +
        labs(x="ACTB", y="GAPDH") +
        theme(axis.text=element_text(size=13), axis.title=element_text(size=15),
            strip.background = element_blank(), strip.text.x = element_blank())
  print(p)
  dev.off()
}
