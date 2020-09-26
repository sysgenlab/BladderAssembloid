### Fig2 c ###
options(java.parameters="-Xmx32g")
library(Seurat)
library(readr)
library(ggplot2)
library(xlsx)
library(pheatmap)
library(cowplot)
library(gridExtra)

Seu.WT <- read_rds("data/scRNA/Seurat/RData/WT.Seu_mm.SCTransform.SNN.Rds")
Seu.Org_day1 <- read_rds("data/scRNA/Seurat/RData/Org_day1.Seu_mm.SCTransform.SNN.Rds")
Seu.Org_day7 <- read_rds("data/scRNA/Seurat/RData/Org_day7.Seu_mm.SCTransform.SNN.Rds")

expr_WT <- as.data.frame(GetAssayData(Seu.WT, slot = "data"))
expr_day7 <- as.data.frame(GetAssayData(Seu.Org_day7, slot = "data"))
expr_day1 <- as.data.frame(GetAssayData(Seu.Org_day1, slot = "data"))

# Cell Types to Show
Type <- c("Epithelial cells", "Stromal fibroblasts", "Smooth muscle cells", "Endothelial cells", "Schwann cells", "Macrophages")

# Function to Extract Expression Matrix
Mat.ExPc <- function(g.list, Type){
  WT <- expr_WT[g.list,Idents(Seu.WT) == Type]
  day1 <- expr_day1[g.list,Idents(Seu.Org_day1) == Type]
  day7 <- expr_day7[g.list,Idents(Seu.Org_day7) == Type]
  
  rownames(WT) <- rownames(day1) <- rownames(day7) <- g.list
  WT[is.na(WT)] <- 0
  day1[is.na(day1)] <- 0
  day7[is.na(day7)] <- 0
  
  expr <- data.frame(WT=apply(WT, 1, function(x) mean(x[x != 0])),
                     day1=apply(day1, 1, function(x) mean(x[x != 0])),
                     day7=apply(day7, 1, function(x) mean(x[x != 0])))
  expr[is.na(expr)] <- 0
  
  return(expr)
}

# Top 50 Marker Genes
top_thresh <- 50

# Legend thresholds
legend <- list(Epi  = seq(0.5, 2.5, length.out=100),
               Str  = seq(0.5, 2.5, length.out=100),
               SMC  = seq(0.5, 2, length.out=100),
               Endo = seq(0.5, 2, length.out=100),
               Sch  = seq(0.5, 1.5, length.out=100),
               Mphi = seq(0.5, 1.5, length.out=100))

legend_breaks <- list(Epi  = seq(0.5, 2.5, 0.5),
                      Str  = seq(0.5, 2.5, 0.5),
                      SMC  = seq(0.5, 2, 0.75),
                      Endo = seq(0.5, 2, 0.75),
                      Sch  = seq(0.5, 1.5, 0.5),
                      Mphi = seq(0.5, 1.5, 0.5))

legend_labels <- list(Epi  = c("0.5", "1", "1.5", "2", ">2.5"),
                      Str  = c("0.5", "1", "1.5", "2", ">2.5"),
                      SMC  = c("<0.5", "1.25", ">2"),
                      Endo = c("0.5", "1.25", ">2"),
                      Sch  = c("<0.5", "1", ">1.5"),
                      Mphi = c("<0.5", "1", ">1.5"))

plot_list=list()
for(i in seq(Type)){
  MG <- read.xlsx2("data/scRNA/Seurat/MG_list/Marker.Genes.Wilcox.WT.mm.Type.xlsx", sheetName=Type[i], stringsAsFactors=FALSE)
  MG <- MG[1:min(top_thresh, nrow(MG)),"gene"]

  expr <- Mat.ExPc(MG, Type[i])

  p <- pheatmap(expr[order(-expr$WT),], show_rownames=TRUE, clustering_distance_cols="euclidean", clustering_method="single", cluster_rows=FALSE,
                 main=Type[i], border_color=NA, breaks=legend[[i]], legend_breaks=legend_breaks[[i]], legend_labels=legend_labels[[i]], silent=TRUE)

  plot_list[[i]] <- p[[4]]
}
g <- grid.arrange(arrangeGrob(grobs=plot_list, nrow=1, widths=c(1.1, 1, 1, 1, 1, 1)))
dev.off()

# The column order should be exchanged manually
pdf("Figure2/image/Fig2c.pdf", width=15, height=6)
plot_grid(g)
dev.off()