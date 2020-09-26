### Fig2 b ###
options(java.parameters="-Xmx32g")
library(Seurat)
library(readr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(RColorBrewer)

Seu.WT <- read_rds("data/scRNA/Seurat/RData/WT.Seu_mm.SCTransform.SNN.Rds")
Seu.Org_day1 <- read_rds("data/scRNA/Seurat/RData/Org_day1.Seu_mm.SCTransform.SNN.Rds")
Seu.Org_day1_hs <- read_rds("data/scRNA/Seurat/RData/Org_day1.Seu_hs.SCTransform.SNN.Rds")
Seu.Org_day7 <- read_rds("data/scRNA/Seurat/RData/Org_day7.Seu_mm.SCTransform.SNN.Rds")
Seu.Org_day7_hs <- read_rds("data/scRNA/Seurat/RData/Org_day7.Seu_hs.SCTransform.SNN.Rds")

levels(Idents(Seu.WT))[levels(Idents(Seu.WT)) == "Endothelial cells"] <- "Endothelial cells (mouse)"
levels(Idents(Seu.Org_day1))[levels(Idents(Seu.Org_day1)) == "Endothelial cells"] <- "Endothelial cells (mouse)"
levels(Idents(Seu.Org_day7))[levels(Idents(Seu.Org_day7)) == "Endothelial cells"] <- "Endothelial cells (mouse)"

df_WT <- data.frame(type=names(summary(Idents(Seu.WT))), count=summary(Idents(Seu.WT)))
df_WT$Proportion <- df_WT$count/sum(df_WT$count)*100
df_WT$Sample <- "WT"

df_Org1 <- data.frame(type=names(summary(Idents(Seu.Org_day1))), count=summary(Idents(Seu.Org_day1)))
df_Org1 <- rbind(df_Org1, data.frame(type="Endothelial cells (human)", count=ncol(Seu.Org_day1_hs)))
df_Org1$Proportion <- df_Org1$count/sum(df_Org1$count)*100
df_Org1$Sample <- "1-day"

df_Org7 <- data.frame(type=names(summary(Idents(Seu.Org_day7))), count=summary(Idents(Seu.Org_day7)))
df_Org7 <- rbind(df_Org7, data.frame(type="Endothelial cells (human)", count=ncol(Seu.Org_day7_hs)))
df_Org7$Proportion <- df_Org7$count/sum(df_Org7$count)*100
df_Org7$Sample <- "7-day"

df <- rbind(df_WT, df_Org1, df_Org7)
lvl <- c("Others (cycling cells)",  
         "Mesothelial cells", 
         "Macrophages", 
         "Schwann cells", 
         "Endothelial cells (human)",
         "Endothelial cells (mouse)", 
         "Smooth muscle cells", 
         "Stromal fibroblasts", 
         "Epithelial cells")
df$type <- factor(df$type, levels=lvl)
df$Sample <- factor(df$Sample, levels=c("WT", "1-day", "7-day"))

p <- ggplot(data=df, aes(x=Sample, y=Proportion, fill=type)) + geom_bar(stat="identity", width=0.7)
colors <- c("#808080",  # Others (cycling cells)
            "#A65628",  # Mesothelial cells
            "#FFFF33",  # Macrophages
            "#FF7F00",  # Schwann cells
            "#CD7DD9",  # Endothelial cells (human)
            "#984EA3",  # Endothelial cells (mouse)
            "#4DAF4A",  # Smooth muscle cells
            "#377EB8",  # Stromal fibroblasts
            "#E41A1C")  # Epithelial cells
p <- p + scale_fill_manual(values=colors)
p <- p + theme_classic()
p <- p + scale_y_continuous(expand=c(0,0))
p <- p + labs(x=NULL, y="Proportion of cells (%)")
p <- p + theme(axis.text.x=element_text(size=11))
p <- p + guides(fill=guide_legend(title=NULL, reverse = TRUE))

pdf("Figure2/image/Fig2b.pdf", width=5, height=5)
print(p)
dev.off()