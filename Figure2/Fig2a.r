### Fig2 a ###
# Epithelial cells: #E41A1C #E61983
# Stromal fibroblasts: #377EB8 #3794B9
# Smooth muscle cells: #4DAF4A #80B04A
# Endothelial cells: #984EA3 #6D4EA2
# Schwann cells: #FF7F00
# Macrophages: #FFFF33
# Mesothelial cells: #A65628
# Others (cycling cells): #808080
options(java.parameters="-Xmx32g")
library(Seurat)
library(readr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(RColorBrewer)

Seu.WT <- read_rds("data/scRNA/Seurat/RData/WT.Seu_mm.SCTransform.SNN.Rds")
Seu.Org_day1 <- read_rds("data/scRNA/Seurat/RData/Org_day1.Seu_mm.SCTransform.SNN.Rds")
Seu.Org_day7 <- read_rds("data/scRNA/Seurat/RData/Org_day7.Seu_mm.SCTransform.SNN.Rds")

pt.size=0.5

p1 <- DimPlot(Seu.WT, reduction="tsne", group.by="Type", pt.size=pt.size) + 
  scale_color_manual(values=c("#E41A1C", 
                              "#377EB8", 
							  "#4DAF4A",
							  "#984EA3", 
							  "#FF7F00", 
                              "#FFFF33", 
                              "#A65628")) + 
  ggtitle("Wild-type bladder") + theme(plot.title=element_text(hjust=0.5))

p2 <- DimPlot(Seu.Org_day1, reduction="tsne", group.by="Type", pt.size=pt.size) + 
  scale_color_manual(values=c("#E41A1C",
                              "#377EB8",
							  "#4DAF4A",
							  "#984EA3",
							  "#FF7F00", 
                              "#FFFF33", 
                              "#808080")) + 
  ggtitle("1-days bladder assembloid") + theme(plot.title=element_text(hjust=0.5))

p3 <- DimPlot(Seu.Org_day7, reduction="tsne", group.by="Type", pt.size=pt.size) + 
  scale_color_manual(values=c("#E41A1C",
                              "#377EB8",
							  "#4DAF4A",
							  "#984EA3", 
							  "#FF7F00", 
                              "#FFFF33")) + 
  ggtitle("7-days bladder assembloid") + theme(plot.title=element_text(hjust=0.5))

p1_legend <- get_legend(p1)
p2_legend <- get_legend(p2)
p3_legend <- get_legend(p3)

p <- plot_grid(p1 + theme(legend.position="none"), p2 + theme(legend.position="none"), p3 + theme(legend.position="none"), nrow=1)
p_legend <- plot_grid(p1_legend, p2_legend, p3_legend, nrow=1)

pdf("Figure2/image/Fig2a.pdf", width=24, height=8)
print(p)
dev.off()

pdf("Figure2/image/Fig2a_legend.pdf", width=24, height=8)
print(p_legend)
dev.off()