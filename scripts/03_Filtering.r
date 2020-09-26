# The script follows previous study: http://oshlacklab.com/combes-organoid-paper/

### Quality Control ###
# Read in single-cell RNA-Seq data, produce various quality control plots and remove any low-quality cells or uniformative genes.
## The number of counts per barcode
## The number of genes per barcode
## The fraction of counts from mitochondrial genes per barcode

library(readr)
library(scater)
library(cowplot)

dir.create("data/scRNA/scater/QCplot_filt")
dir.create("data/scRNA/scater/RData_filt")

sample <- c("WT", "Org_day1", "Org_day7")

for(i in 1:length(sample)){
  sce_mm <- read_rds(paste0("data/scRNA/scater/RData/", sample[i], ".sce_mm.QC.norm.Rds"))
  sce_hs <- read_rds(paste0("data/scRNA/scater/RData/", sample[i], ".sce_hs.QC.norm.Rds"))

  colData(sce_mm)$Filtered <- FALSE
  colData(sce_hs)$Filtered <- FALSE

  dir.create(paste0("data/scRNA/scater/QCplot_filt/", sample[i]))

  # Total Counts: Depth
  ##scater-1.12.2 - plotColData
  thresh.count <- 60000
  
  pdf(paste0("data/scRNA/scater/QCplot_filt/", sample[i], "/01_Total_Counts.mm.pdf"), height=7, width=5)
  p <- plotColData(sce_mm, x = "Sample", y = "total_counts") +
        labs(x="Human", y="Count depth") + scale_x_discrete(labels = sample[i]) +
        theme(axis.text=element_text(size=13), axis.title=element_text(size=15)) +
        geom_hline(yintercept=thresh.count, colour = "red", size=1.5, linetype="dashed")
  print(p)
  dev.off()
 
  pdf(paste0("data/scRNA/scater/QCplot_filt/", sample[i], "/01_Total_Counts.hs.pdf"), height=7, width=5)
  p <- plotColData(sce_hs, x = "Sample", y = "total_counts") +
        labs(x="Human", y="Count depth") + scale_x_discrete(labels = sample[i]) +
        theme(axis.text=element_text(size=13), axis.title=element_text(size=15)) +
        geom_hline(yintercept=thresh.count, colour = "red", size=1.5, linetype="dashed")
  print(p)
  dev.off()
 
  colData(sce_mm)$Filtered <- colData(sce_mm)$Filtered | colData(sce_mm)$total_counts > thresh.count
  colData(sce_hs)$Filtered <- colData(sce_hs)$Filtered | colData(sce_hs)$total_counts > thresh.count
    
  # Number of Features by Library Size
  ##scater-1.12.2 - plotColData
  ## Cells that express many genes are potential multiplets
  thresh.gene <- c(500, 8000)
  
  pdf(paste0("data/scRNA/scater/QCplot_filt/", sample[i], "/02_Features_by_Library_Size.mm.pdf"), height=7, width=9)
  p <- plotColData(sce_mm, x = "total_counts", y = "total_features_by_counts", colour_by = "pct_dropout", shape_by="Filtered") +
        labs(x="Count depth", y="Number of genes") +
        theme(axis.text=element_text(size=13), axis.title=element_text(size=15)) +
        geom_vline(xintercept=thresh.count, colour = "red", size=1.5, linetype="dashed") + 
        geom_hline(yintercept=thresh.gene, colour = "red", size=1.5, linetype="dashed")
  print(p)
  dev.off() 

  pdf(paste0("data/scRNA/scater/QCplot_filt/", sample[i], "/02_Features_by_Library_Size.hs.pdf"), height=7, width=9)
  p <- plotColData(sce_hs, x = "total_counts", y = "total_features_by_counts", colour_by = "pct_dropout", shape_by="Filtered") +
        labs(x="Count depth", y="Number of genes") +
        theme(axis.text=element_text(size=13), axis.title=element_text(size=15)) +
        geom_vline(xintercept=thresh.count, colour = "red", size=1.5, linetype="dashed") +
        geom_hline(yintercept=thresh.gene, colour = "red", size=1.5, linetype="dashed")
  print(p)
  dev.off()
 
  colData(sce_mm)$Filtered <- colData(sce_mm)$Filtered | colData(sce_mm)$total_features_by_counts < thresh.gene[1]
  colData(sce_mm)$Filtered <- colData(sce_mm)$Filtered | colData(sce_mm)$total_features_by_counts > thresh.gene[2]
  colData(sce_hs)$Filtered <- colData(sce_hs)$Filtered | colData(sce_hs)$total_features_by_counts < thresh.gene[1]
  colData(sce_hs)$Filtered <- colData(sce_hs)$Filtered | colData(sce_hs)$total_features_by_counts > thresh.gene[2]
    
    
  # Mitochondrial genes
  ## Overexpression of mitochondrial genes can be an indication that cell is stressed or damaged in some way.
  thresh.MT <- 30
  
  pdf(paste0("data/scRNA/scater/QCplot_filt/", sample[i], "/03_Mitochondrial_Genes.mm.pdf"), height=7, width=5)
  p <- plotColData(sce_mm, x = "Sample", y = "PctCountsMT", shape_by="Filtered") +
        labs(x="Mouse", y="% counts MT") + scale_x_discrete(labels = sample[i]) +
        theme(axis.text=element_text(size=13), axis.title=element_text(size=15)) +
        geom_hline(yintercept=thresh.MT, colour = "red", size=1.5, linetype="dashed")
  print(p)
  dev.off()

  pdf(paste0("data/scRNA/scater/QCplot_filt/", sample[i], "/03_Mitochondrial_Genes.hs.pdf"), height=7, width=5)
  p <- plotColData(sce_hs, x = "Sample", y = "PctCountsMT", shape_by="Filtered") +
        labs(x="Human", y="% counts MT") + scale_x_discrete(labels = sample[i]) +
        theme(axis.text=element_text(size=13), axis.title=element_text(size=15)) +
        geom_hline(yintercept=thresh.MT, colour = "red", size=1.5, linetype="dashed")
  print(p)
  dev.off()
 
  colData(sce_mm)$Filtered <- colData(sce_mm)$Filtered | colData(sce_mm)$PctCountsMT > thresh.MT
  colData(sce_hs)$Filtered <- colData(sce_hs)$Filtered | colData(sce_hs)$PctCountsMT > thresh.MT
  
  
  # Ribosomal genes
  thresh.ribo <- 45
  
  pdf(paste0("data/scRNA/scater/QCplot_filt/", sample[i], "/04_Ribosomal_Genes.mm.pdf"), height=7, width=5)
  p <- plotColData(sce_mm, x = "Sample", y = "PctCountsRibo", shape_by="Filtered") +
        labs(x="Mouse", y="% counts risomal") + scale_x_discrete(labels = sample[i]) +
        theme(axis.text=element_text(size=13), axis.title=element_text(size=15)) +
        geom_hline(yintercept=thresh.ribo, colour = "red", size=1.5, linetype="dashed")
  print(p)
  dev.off()

  pdf(paste0("data/scRNA/scater/QCplot_filt/", sample[i], "/04_Ribosomal_Genes.hs.pdf"), height=7, width=5)
  p <- plotColData(sce_hs, x = "Sample", y = "PctCountsRibo", shape_by="Filtered") +
        labs(x="Human", y="% counts risomal") + scale_x_discrete(labels = sample[i]) +
        theme(axis.text=element_text(size=13), axis.title=element_text(size=15)) +
        geom_hline(yintercept=thresh.ribo, colour = "red", size=1.5, linetype="dashed")
  print(p)
  dev.off()
  
  colData(sce_mm)$Filtered <- colData(sce_mm)$Filtered | colData(sce_mm)$PctCountsRibo > thresh.ribo
  colData(sce_hs)$Filtered <- colData(sce_hs)$Filtered | colData(sce_hs)$PctCountsRibo > thresh.ribo
    
    
  # Housekeeping genes
  thresh.gapdh <- 0.5
  thresh.actb <- 0.5
  
  actb.id <- rowData(sce_mm)[which(rowData(sce_mm)$symbol == "Actb"),1]
  gapdh.id <- rowData(sce_mm)[which(rowData(sce_mm)$symbol == "Gapdh"),1]
  key <- c("Actb", "Gapdh")
  names(key) <- c(actb.id, gapdh.id)

  pdf(paste0("data/scRNA/scater/QCplot_filt/", sample[i], "/05_Housekeeping_Genes.mm.pdf"), height=7, width=8)
  p <- plotExpression(sce_mm, gapdh.id, x = actb.id, shape_by="Filtered") +
        labs(x="Actb", y="Gapdh") +
        geom_hline(yintercept = thresh.gapdh, colour = "red", size = 1.5, linetype = "dashed") +
        geom_vline(xintercept = thresh.actb, colour = "red", size = 1.5, linetype = "dashed") +
        theme(axis.text=element_text(size=13), axis.title=element_text(size=15),
            strip.background = element_blank(), strip.text.x = element_blank())
  print(p)
  dev.off()

  colData(sce_mm)$Filtered <- colData(sce_mm)$Filtered | exprs(sce_mm)[actb.id, ] < thresh.actb | exprs(sce_mm)[gapdh.id, ] < thresh.gapdh

  actb.id <- rowData(sce_hs)[which(rowData(sce_hs)$symbol == "ACTB"),1]
  gapdh.id <- rowData(sce_hs)[which(rowData(sce_hs)$symbol == "GAPDH"),1]
  key <- c("ACTB", "GAPDH")
  names(key) <- c(actb.id, gapdh.id)

  pdf(paste0("data/scRNA/scater/QCplot_filt/", sample[i], "/05_Housekeeping_Genes.hs.pdf"), height=7, width=8)
  p <- plotExpression(sce_hs, gapdh.id, x = actb.id, shape_by="Filtered") +
        labs(x="ACTB", y="GAPDH") +
        geom_hline(yintercept = thresh.gapdh, colour = "red", size = 1.5, linetype = "dashed") +
        geom_vline(xintercept = thresh.actb, colour = "red", size = 1.5, linetype = "dashed") +
        theme(axis.text=element_text(size=13), axis.title=element_text(size=15),
            strip.background = element_blank(), strip.text.x = element_blank())
  print(p)
  dev.off()

  colData(sce_hs)$Filtered <- colData(sce_hs)$Filtered | exprs(sce_hs)[actb.id, ] < thresh.actb | exprs(sce_hs)[gapdh.id, ] < thresh.gapdh
  
  
  # Filtering Barcodes
  sce_mm <- sce_mm[, !colData(sce_mm)$Filtered]
  sce_hs <- sce_hs[, !colData(sce_hs)$Filtered]
  
  
  # Gene Expression
  ## Genes that have less than two counts across all cells (not expressed genes)
  keep <- rowSums(counts(sce_mm)) > 1
  sce_mm <- sce_mm[keep,]
  keep <- rowSums(counts(sce_hs)) > 1
  sce_hs <- sce_hs[keep,]
  
  ## Genes that are expressed in less than two individual cells
  keep <- rowSums(counts(sce_mm) != 0) > 1
  sce_mm <- sce_mm[keep, ]
  keep <- rowSums(counts(sce_hs) != 0) > 1
  sce_hs <- sce_hs[keep, ]  

  # MGI genes
  ## Genes that don't have MGI symbols mostly pseudogenes and are unlikely to be informative.
  keep <- !(rowData(sce_mm)$mgi_symbol == "") & !(is.na(rowData(sce_mm)$mgi_symbol))
  sce_mm <- sce_mm[keep, ]

  keep <- !(rowData(sce_hs)$hgnc_symbol == "") & !(is.na(rowData(sce_hs)$hgnc_symbol))
  sce_hs <- sce_hs[keep, ]
  
  dups <- which(duplicated(rowData(sce_mm)$symbol))
  rowData(sce_mm)[dups, "symbol"] <- paste0(rowData(sce_mm)[dups, "symbol"], "-dup")

  dups <- which(duplicated(rowData(sce_hs)$symbol))
  rowData(sce_hs)[dups, "symbol"] <- paste0(rowData(sce_hs)[dups, "symbol"], "-dup")
  
  write_rds(sce_mm, paste0("data/scRNA/scater/RData_filt/", sample[i], ".sce_mm.QC.norm.filt.Rds"))
  write_rds(sce_hs, paste0("data/scRNA/scater/RData_filt/", sample[i], ".sce_hs.QC.norm.filt.Rds"))
}
