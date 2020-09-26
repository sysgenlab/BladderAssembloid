### Extended Fig3 a ###
# HULEC: #CD7DD9
options(java.parameters="-Xmx32g")
library(Seurat)
library(readr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(RColorBrewer)

Seu.Org_day1_hs <- read_rds("data/scRNA/Seurat/RData/Org_day1.Seu_hs.SCTransform.SNN.Rds")
Seu.Org_day7_hs <- read_rds("data/scRNA/Seurat/RData/Org_day7.Seu_hs.SCTransform.SNN.Rds")

pt.size=0.5

p1 <- DimPlot(Seu.Org_day1_hs, reduction="tsne", group.by="Type", pt.size=pt.size) + 
  scale_color_manual(values="#CD7DD9") + 
  ggtitle("1-day") + theme(plot.title=element_text(hjust=0.5))

p2 <- DimPlot(Seu.Org_day7_hs, reduction="tsne", group.by="Type", pt.size=pt.size) + 
  scale_color_manual(values="#CD7DD9") + 
  ggtitle("7-days") + theme(plot.title=element_text(hjust=0.5))

p <- plot_grid(p1, p2, ncol=1)

pdf("Extended_Figure3/image/Extended_Fig3a.pdf", width=4, height=8)
print(p)
dev.off()