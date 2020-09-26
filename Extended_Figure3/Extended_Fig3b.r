### Endothelial cells and HULEC: Expanded Fig 3b ###
options(java.parameters="-Xmx32g")
library(Seurat)
library(readr)
library(ggplot2)
library(xlsx)
library(RColorBrewer)

Seu.WT <- read_rds("data/scRNA/Seurat/RData/WT.Seu_mm.SCTransform.SNN.Rds")
Seu.Org_day1_hs <- read_rds("data/scRNA/Seurat/RData/Org_day1.Seu_hs.SCTransform.SNN.Rds")
Seu.Org_day7_hs <- read_rds("data/scRNA/Seurat/RData/Org_day7.Seu_hs.SCTransform.SNN.Rds")

mart <- read.delim("data/mart_export.txt")
mart_ortho <- unique(mart[,c("Gene.name", "Mouse.gene.name")])
mart_ortho <- mart_ortho[mart_ortho$Gene.name != "" & mart_ortho$Mouse.gene.name != "",]

MG.1 <- as.character(read.xlsx2("data/scRNA/Seurat/MG_list/Marker.Genes.Wilcox.WT.mm.xlsx", sheetName="cluster3")$gene)
MG.2 <- as.character(read.xlsx2("data/scRNA/Seurat/MG_list/Marker.Genes.Wilcox.WT.mm.xlsx", sheetName="cluster5")$gene)

# Org day1 HULEC
expr_Org1_hs <- as.matrix(GetAssayData(Seu.Org_day1_hs, slot = "data"))
ortho_Org1_hs <- as.character(mart_ortho$Mouse.gene.name[match(rownames(expr_Org1_hs), mart_ortho$Gene.name)])

# Org day7 HULEC
expr_Org7_hs <- as.matrix(GetAssayData(Seu.Org_day7_hs, slot = "data"))
ortho_Org7_hs <- as.character(mart_ortho$Mouse.gene.name[match(rownames(expr_Org7_hs), mart_ortho$Gene.name)])

# WT Endothelial cells 1
expr_WT_1 <- as.matrix(GetAssayData(Seu.WT[,Seu.WT@meta.data$Subtype == "Endothelial cells 1"], slot = "data"))

# WT Endothelial cells 2
expr_WT_2 <- as.matrix(GetAssayData(Seu.WT[,Seu.WT@meta.data$Subtype == "Endothelial cells 2"], slot = "data"))


df_1 <- data.frame()
MG <- character(0)
i <- 1
while(length(MG) < 15){
  gene <- MG.1[i]
  if(sum(ortho_Org1_hs == gene, na.rm=TRUE) != 0 & sum(ortho_Org7_hs == gene, na.rm=TRUE) != 0){
    ortholog <- as.character(mart_ortho$Gene.name[mart_ortho$Mouse.gene.name == gene])[1]
	
	WT_idx <- expr_WT_1[gene,] != 0
	Org1_idx <- expr_Org1_hs[ortholog,] != 0
	Org7_idx <- expr_Org7_hs[ortholog,] != 0
	
	pct_WT <- sum(WT_idx)/ncol(expr_WT_1)*100
	pct_Org1 <- sum(Org1_idx)/ncol(expr_Org1_hs)*100
	pct_Org7 <- sum(Org7_idx)/ncol(expr_Org7_hs)*100
	
	avg_WT <- mean(expr_WT_1[gene, WT_idx])
	avg_Org1 <- mean(expr_Org1_hs[ortholog, Org1_idx])
	avg_Org7 <- mean(expr_Org7_hs[ortholog, Org7_idx])
    
	df_tmp <- data.frame(x=c("WT\nA3", "1d\nhB0", "7d\nhC0"), y=gene, pct=c(pct_WT, pct_Org1, pct_Org7), avg=c(avg_WT, avg_Org1, avg_Org7))
    df_1 <- rbind(df_1, df_tmp)
  
    MG <- c(gene, MG)
  }
  i <- i+1
}

df_1$x <- factor(df_1$x, levels=c("WT\nA3", "1d\nhB0", "7d\nhC0"))
df_1$y <- factor(df_1$y, levels=MG)
df_1$Type <- "A3"

mart_ortho <- mart_ortho[-8505,]
df_2 <- data.frame()
MG <- character(0)
i <- 1
while(length(MG) < 15){
  gene <- MG.2[i]
  if(sum(ortho_Org1_hs == gene, na.rm=TRUE) != 0 & sum(ortho_Org7_hs == gene, na.rm=TRUE) != 0){
    ortholog <- as.character(mart_ortho$Gene.name[mart_ortho$Mouse.gene.name == gene])[1]
	
	WT_idx <- expr_WT_2[gene,] != 0
	Org1_idx <- expr_Org1_hs[ortholog,] != 0
	Org7_idx <- expr_Org7_hs[ortholog,] != 0
	
	pct_WT <- sum(WT_idx)/ncol(expr_WT_2)*100
	pct_Org1 <- sum(Org1_idx)/ncol(expr_Org1_hs)*100
	pct_Org7 <- sum(Org7_idx)/ncol(expr_Org7_hs)*100
	
	avg_WT <- mean(expr_WT_2[gene, WT_idx])
	avg_Org1 <- mean(expr_Org1_hs[ortholog, Org1_idx])
	avg_Org7 <- mean(expr_Org7_hs[ortholog, Org7_idx])
    
	df_tmp <- data.frame(x=c("WT\nA5", "1d\nhB0", "7d\nhC0"), y=gene, pct=c(pct_WT, pct_Org1, pct_Org7), avg=c(avg_WT, avg_Org1, avg_Org7))
    df_2 <- rbind(df_2, df_tmp)
  
    MG <- c(gene, MG)
  }
  i <- i+1
}

df_2$x <- factor(df_2$x, levels=c("WT\nA5", "1d\nhB0", "7d\nhC0"))
df_2$y <- factor(df_2$y, levels=MG)
df_2$Type <- "A5"

df <- rbind(df_1, df_2)
df$Type <- factor(df$Type, levels=c("A3", "A5"))
df$x <- factor(df$x, levels=c("WT\nA3", "WT\nA5", "1d\nhB0", "7d\nhC0"))
levels(df$Type) <- c("Endothelial cells\n(WT A3 & HULEC)", "Endothelial cells\n(WT A5 & HULEC)")

p <- ggplot(df, aes(x=x, y=y)) + geom_tile(aes(fill=avg), colour="black", size=0.25)
p <- p + facet_wrap(~ Type, scales="free", nrow=1)
p <- p + scale_fill_gradient(low="purple", high="gold")
p <- p + theme_bw()
p <- p + labs(x=NULL, y=NULL)
p <- p + theme(axis.title=element_blank(), 
               axis.text=element_text(size=12), 
               strip.text=element_text(size=15), 
               strip.background=element_blank(),
               axis.ticks=element_blank(),
               panel.border = element_rect(size=0.5))
p <- p + scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0))

pdf("Extended_Figure3/image/Extended_Fig3b.pdf", width=5, height=7.5)
print(p)
dev.off()