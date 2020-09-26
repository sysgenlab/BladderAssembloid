### Subtype expression comparison: Expanded Fig 3d ###
options(java.parameters="-Xmx32g")
library(Seurat)
library(readr)
library(ggplot2)
library(xlsx)
library(RColorBrewer)

Seu.WT <- read_rds("data/scRNA/Seurat/RData/WT.Seu_mm.SCTransform.SNN.Rds")
Seu.Org_day1 <- read_rds("data/scRNA/Seurat/RData/Org_day1.Seu_mm.SCTransform.SNN.Rds")
Seu.Org_day7 <- read_rds("data/scRNA/Seurat/RData/Org_day7.Seu_mm.SCTransform.SNN.Rds")

smpl <- c("Epithelial cells", "Stromal fibroblasts", "Smooth muscle cells", "Endothelial cells")

dff <- data.frame()
for(i in 1:4){
  MG <- as.character(read.xlsx2("data/scRNA/Seurat/MG_list/Marker.Genes.Wilcox.WT.mm.Type.xlsx", sheetName=smpl[i])$gene)

  A <- unique(Seu.WT@meta.data[Seu.WT@meta.data$Type == smpl[i],"seurat_clusters"])
  B <- unique(Seu.Org_day1@meta.data[Seu.Org_day1@meta.data$Type == smpl[i],"seurat_clusters"])
  C_1 <- unique(Seu.Org_day7@meta.data[Seu.Org_day7@meta.data$Subtype == paste0(smpl[i], " 1"),"seurat_clusters"])
  C_2 <- unique(Seu.Org_day7@meta.data[Seu.Org_day7@meta.data$Subtype == paste0(smpl[i], " 2"),"seurat_clusters"])
  
  A <- A[order(A)][1]
  B <- B[order(B)][1]

  expr_WT <- as.matrix(GetAssayData(Seu.WT[,Seu.WT@meta.data$seurat_clusters == A], slot = "data"))
  expr_Org1 <- as.matrix(GetAssayData(Seu.Org_day1[,Seu.Org_day1@meta.data$seurat_clusters == B], slot = "data"))
  expr_Org7_1 <- as.matrix(GetAssayData(Seu.Org_day7[,Seu.Org_day7@meta.data$seurat_clusters == C_1], slot = "data"))
  expr_Org7_2 <- as.matrix(GetAssayData(Seu.Org_day7[,Seu.Org_day7@meta.data$seurat_clusters == C_2], slot = "data"))
  
  df <- data.frame()
  for(j in 1:15){
    gene <- MG[j]
  
    idx_expr_WT <- expr_WT[gene,] != 0
    idx_expr_Org1 <- expr_Org1[gene,] != 0
    idx_expr_Org7_1 <- expr_Org7_1[gene,] != 0
	idx_expr_Org7_2 <- expr_Org7_2[gene,] != 0
  
    pct_WT <- sum(idx_expr_WT)/ncol(expr_WT)*100
    pct_Org7_1 <- sum(idx_expr_Org7_1)/ncol(expr_Org7_1)*100
    pct_Org7_2 <- sum(idx_expr_Org7_2)/ncol(expr_Org7_2)*100
	pct_Org1 <- sum(idx_expr_Org1)/ncol(expr_Org1)*100
  
    avg_WT <- mean(expr_WT[gene,idx_expr_WT])
    avg_Org7_1 <- mean(expr_Org7_1[gene,idx_expr_Org7_1])
    avg_Org7_2 <- mean(expr_Org7_2[gene,idx_expr_Org7_2])
	avg_Org1 <- mean(expr_Org1[gene,idx_expr_Org1])
  
    df_tmp <- data.frame(x=c(paste0("WT\nA", A), paste0("1d\nB", B), paste0("7d\nC", C_1), paste0("7d\nC", C_2)), 
                         y=gene, 
                         pct=c(pct_WT, pct_Org1, pct_Org7_1, pct_Org7_2), 
                         avg=c(avg_WT, avg_Org1, avg_Org7_1, avg_Org7_2))
    df <- rbind(df, df_tmp)
  }

  df$y <- factor(df$y, levels=MG[15:1])
  
  df$Type <- smpl[i]
  dff <- rbind(dff, df)
}

Type <- c("Epithelial\ncells", "Stromal\nfibroblasts", "Smooth\nmuscle cells", "Endothelial\ncells")
names(Type) <- smpl[1:4]
dff$Type <- factor(Type[dff$Type], levels=Type)
dff$x <- relevel(dff$x, "7d\nC3")
dff$x <- relevel(dff$x, "1d\nB2")

p <- ggplot(dff, aes(x=x, y=y)) + geom_tile(aes(fill=avg), colour="black", size=0.25)
p <- p + facet_wrap(~ Type, scales="free", nrow=1)
p <- p + scale_fill_gradient(low="purple", high="gold")
p <- p + theme_bw()
p <- p + labs(x=NULL, y=NULL)
p <- p + theme(axis.title=element_blank(), 
               axis.text=element_text(size=10), 
			   strip.text=element_text(size=12), 
			   strip.background=element_blank(), 
			   axis.ticks=element_blank(),
			   panel.border = element_rect(size=0.5))
p <- p + scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0))

pdf("Extended_Figure3/image/Extended_Fig3d.pdf", width=7, height=4.5)
print(p)
dev.off()