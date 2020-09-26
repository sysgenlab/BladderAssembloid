### Endothelial cell subtype comparison: Expanded Fig 3c ###
options(java.parameters="-Xmx32g")
library(Seurat)
library(readr)
library(ggplot2)
library(cowplot)
library(xlsx)
library(RColorBrewer)

Seu.WT <- read_rds("data/scRNA/Seurat/RData/WT.Seu_mm.SCTransform.SNN.Rds")
Seu.Org_day1 <- read_rds("data/scRNA/Seurat/RData/Org_day1.Seu_mm.SCTransform.SNN.Rds")
Seu.Org_day7 <- read_rds("data/scRNA/Seurat/RData/Org_day7.Seu_mm.SCTransform.SNN.Rds")

expr_Org1 <- as.matrix(GetAssayData(Seu.Org_day1[,Seu.Org_day1@meta.data$Subtype == "Endothelial cells"], slot = "data"))
expr_Org7_1 <- as.matrix(GetAssayData(Seu.Org_day7[,Seu.Org_day7@meta.data$Subtype == "Endothelial cells 1"], slot = "data"))
expr_Org7_2 <- as.matrix(GetAssayData(Seu.Org_day7[,Seu.Org_day7@meta.data$Subtype == "Endothelial cells 2"], slot = "data"))

smpl <- c("Endothelial cells 1", "Endothelial cells 2")
names(smpl) <- c("cluster3", "cluster5")
WT_sub <- c("A3", "A5")

for(i in seq(smpl)){
  MG <- as.character(read.xlsx2("data/scRNA/Seurat/MG_list/Marker.Genes.Wilcox.WT.mm.xlsx", sheetName=names(smpl)[i])$gene)

  expr_WT <- as.matrix(GetAssayData(Seu.WT[,Seu.WT@meta.data$Subtype == smpl[i]], slot = "data"))

  df <- data.frame()
  for(j in 1:15){
    gene <- MG[j]
  
    if(gene %in% rownames(expr_Org1)){
      idx_expr_WT <- expr_WT[gene,] != 0
      idx_expr_Org1 <- expr_Org1[gene,] != 0
      idx_expr_Org7_1 <- expr_Org7_1[gene,] != 0
	  idx_expr_Org7_2 <- expr_Org7_2[gene,] != 0
  
      pct_WT <- sum(idx_expr_WT)/ncol(expr_WT)*100
      pct_Org1 <- sum(idx_expr_Org1)/ncol(expr_Org1)*100
      pct_Org7_1 <- sum(idx_expr_Org7_1)/ncol(expr_Org7_1)*100
	  pct_Org7_2 <- sum(idx_expr_Org7_2)/ncol(expr_Org7_2)*100
  
      avg_WT <- mean(expr_WT[gene,idx_expr_WT])
      avg_Org1 <- mean(expr_Org1[gene,idx_expr_Org1])
      avg_Org7_1 <- mean(expr_Org7_1[gene,idx_expr_Org7_1])
	  avg_Org7_2 <- mean(expr_Org7_2[gene,idx_expr_Org7_2])
	}else{
	  pct_Org1 <- avg_Org1 <- 0
	  
	  idx_expr_WT <- expr_WT[gene,] != 0
      idx_expr_Org7_1 <- expr_Org7_1[gene,] != 0
	  idx_expr_Org7_2 <- expr_Org7_2[gene,] != 0
	
	  pct_WT <- sum(idx_expr_WT)/ncol(expr_WT)*100
      pct_Org7_1 <- sum(idx_expr_Org7_1)/ncol(expr_Org7_1)*100
	  pct_Org7_2 <- sum(idx_expr_Org7_2)/ncol(expr_Org7_2)*100
	  
	  avg_WT <- mean(expr_WT[gene,idx_expr_WT])
      avg_Org7_1 <- mean(expr_Org7_1[gene,idx_expr_Org7_1])
	  avg_Org7_2 <- mean(expr_Org7_2[gene,idx_expr_Org7_2])
	}
  
    df_tmp <- data.frame(x=c(paste0("WT\n", WT_sub[i]), "1d\nB4", "7d\nC6", "7d\nC7"), 
                         y=gene, pct=c(pct_WT, pct_Org1, pct_Org7_1, pct_Org7_2), 
                         avg=c(avg_WT, avg_Org1, avg_Org7_1, avg_Org7_2))
    df <- rbind(df, df_tmp)
  }

  df$x <- factor(df$x, levels=c(paste0("WT\n", WT_sub[i]), "1d\nB4", "7d\nC6", "7d\nC7"))
  df$y <- factor(df$y, levels=MG[15:1])

  assign(paste0("df_", i), df)
}

size <- 6
color <- brewer.pal(3,"Set1")[c(1,3,2,2)]

p1 <- ggplot(df_1, aes(x=x, y=y)) + geom_point(aes(size=pct, color=x, alpha=avg))
p1 <- p1 + scale_alpha_continuous(range=c(0.3,1)) + scale_size_continuous(range=c(1,size))
p1 <- p1 + scale_color_manual(values=color)
p1 <- p1 + theme_classic(base_line_size=1)
p1 <- p1 + labs(x=NULL, y=NULL)
p1 <- p1 + theme(axis.title=element_blank(), 
                axis.text=element_text(size=10), 
                strip.text=element_text(size=12), 
                strip.background=element_blank(),
                legend.position="none")

df_2$y <- relevel(df_2$y, "Cldn5")
p2 <- ggplot(df_2, aes(x=x, y=y)) + geom_point(aes(size=pct, color=x, alpha=avg))
p2 <- p2 + scale_alpha_continuous(range=c(0.3,1)) + scale_size_continuous(range=c(1,size))
p2 <- p2 + scale_color_manual(values=color)
p2 <- p2 + theme_classic(base_line_size=1)
p2 <- p2 + labs(x=NULL, y=NULL)
p2 <- p2 + theme(axis.title=element_blank(), 
                axis.text=element_text(size=10), 
                strip.text=element_text(size=12), 
                strip.background=element_blank(),
                legend.position="none")

main <- "Comparison of endothelial cells\nwith WT A3 (left) and WT A5 (right)"
title <- ggdraw() + draw_label(main, x=0, hjust=0.5) + theme(plot.margin = margin(0,0,0,130))
p <- plot_grid(title, plot_grid(p1, p2), ncol=1, rel_heights=c(0.1,1))

pdf("Extended_Figure3/image/Extended_Fig3c.pdf", width=3.5, height=5)
print(p)
dev.off()