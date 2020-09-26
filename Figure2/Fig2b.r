### Fig2 d ###
options(java.parameters="-Xmx32g")
library(Seurat)
library(readr)
library(ggplot2)
library(xlsx)
library(pheatmap)
library(cowplot)
library(gridExtra)
library(grid)

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

# Marker Gene p-value threshold
pval <- 0.05

# Setting Minimum Correlation
min <- 1  
for(i in seq(Type)){
  MG <- read.xlsx2("data/scRNA/Seurat/MG_list/Marker.Genes.Wilcox.WT.mm.Type.xlsx", sheetName=Type[i], stringsAsFactors=FALSE)
  MG <- MG[as.numeric(MG$p_val_adj) < pval,"gene"]
  
  expr <- Mat.ExPc(MG, Type[i])
  
  res <- cor(expr)
  min <- min(min, min(res))
}

plot_list=list()
for(i in seq(Type)){
  MG <- read.xlsx2("data/scRNA/Seurat/MG_list/Marker.Genes.Wilcox.WT.mm.Type.xlsx", sheetName=Type[i], stringsAsFactors=FALSE)
  MG <- MG[as.numeric(MG$p_val_adj) < pval,"gene"]
    
  expr <- Mat.ExPc(MG, Type[i])
  
  res <- cor(expr)
  mat <- t(res[1,])
  rownames(mat) <- "WT"
  colnames(mat) <- c("WT", "1-day", "7-day")
  p <- pheatmap(mat, display_numbers=TRUE, number_color="black", border_color=NA, cluster_cols=FALSE, cluster_rows=FALSE, 
                main=Type[i], fontsize_number=13, fontsize=13, number_format = "%.3f", breaks=seq(min, 1, length.out=100),
                angle_col=0, legend=FALSE, silent=TRUE)
    
  plot_list[[i]] <- p[[4]]
}
g <- grid.arrange(arrangeGrob(grobs=plot_list, ncol=2))

p <- pheatmap(mat, display_numbers=TRUE, number_color="black", border_color=NA, cluster_cols=FALSE, cluster_rows=FALSE, 
              main=Type[i], fontsize_number=13, fontsize=13, number_format = "%.3f", breaks=seq(min, 1, length.out=100),
              angle_col=0, legend=TRUE, silent=TRUE)

pdf("Figure2/image/Fig2d.pdf", width=7, height=5)
plot_grid(g)
dev.off()

pdf("Figure2/image/Fig2d_legend.pdf", width=7, height=5)
grid.draw(p$gtable$grob[[5]])
dev.off()