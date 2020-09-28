### Extended Fig3 i ###
options(java.parameters="-Xmx256g")
library(Seurat)
library(readr)
library(ggplot2)
library(cowplot)

Seu.Org_day1 <- read_rds("data/scRNA/Seurat/RData/Org_day1.Seu_mm.Epi.SCTransform.SNN.Rds")
Seu.Org_day7 <- read_rds("data/scRNA/Seurat/RData/Org_day7.Seu_mm.Epi.SCTransform.SNN.Rds")

# Rename Clusters
levels(Idents(Seu.Org_day1)) <- c("BP1", "BP2", "B", "I", "L")
levels(Idents(Seu.Org_day7)) <- c("B", "IL", "IH", "L")

genes <- c("Krt5", "Trp63", "Krt8", "Krt18", "Upk1a", "Upk1b", "Upk2", "Upk3a")

color_1 <- c("grey 50", "grey80", "#00BA38", "#619CFF", "#F8766D")
color_2 <- c("#00BA38", "#61D0FF", "#6166FF", "#F8766D")

for(i in seq(genes)){
  p <- VlnPlot(Seu.Org_day1, features=genes[i], pt.size=0) 
  p <- p + scale_fill_manual(values=color_1) 
  p <- p + labs(x=NULL)
  p <- p + theme(legend.position="none", 
                 axis.text.x=element_text(angle=0, hjust=0.5))

  assign(paste0("p", i), p)
}

title <- ggdraw() + draw_label("1-day bladder assembloid", fontface="bold", size=17, hjust=0.5)
q1 <- plot_grid(title, plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, ncol=4), ncol=1, rel_heights=c(0.1,1))

for(i in seq(genes)){
  p <- VlnPlot(Seu.Org_day7, features=genes[i], pt.size=0) 
  p <- p + scale_fill_manual(values=color_2) 
  p <- p + labs(x=NULL)
  p <- p + theme(legend.position="none", 
                 axis.text.x=element_text(angle=0, hjust=0.5))

  assign(paste0("p", i), p)
}

title <- ggdraw() + draw_label("7-days bladder assembloid", fontface="bold", size=17, hjust=0.5)
q2 <- plot_grid(title, plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, ncol=4), ncol=1, rel_heights=c(0.1,1))

q <- plot_grid(q1, q2, ncol=1)

pdf("Extended_Figure3/image/Extended_Fig3i.pdf", width=9, height=10)
print(q)
dev.off()
