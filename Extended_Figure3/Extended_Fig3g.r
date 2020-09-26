### Extended Fig3 e ###
# Basal proliferative 1: grey 50
# Basal proliferative 2: grey 80
# Basal: #00BA38
# Intermediate: #619CFF
# Intermediate low: #61D0FF
# Intermediate high: #6166FF
# Luminal: #F8766D
options(java.parameters="-Xmx256g")
library(Seurat)
library(readr)
library(ggplot2)
library(cowplot)

Seu.Org_day1 <- read_rds("data/scRNA/Seurat/RData/Org_day1.Seu_mm.Epi.SCTransform.SNN.Rds")
Seu.Org_day7 <- read_rds("data/scRNA/Seurat/RData/Org_day7.Seu_mm.Epi.SCTransform.SNN.Rds")

color_1 <- c("grey 50", "grey80", "#00BA38", "#619CFF", "#F8766D")
color_2 <- c("#00BA38", "#61D0FF", "#6166FF", "#F8766D")

pt.size=0.5

p1 <- DimPlot(Seu.Org_day1, reduction="tsne", group.by="Type", pt.size=pt.size) + 
  scale_color_manual(values=color_1) + 
  ggtitle("1-day bladder assembloid") + theme(plot.title=element_text(hjust=0.5))

p2 <- DimPlot(Seu.Org_day7, reduction="tsne", group.by="Type", pt.size=pt.size) + 
  scale_color_manual(values=color_2) + 
  ggtitle("7-days bladder assembloid") + theme(plot.title=element_text(hjust=0.5))

p1_legend <- get_legend(p1)
p2_legend <- get_legend(p2)

p <- plot_grid(p1 + theme(legend.position="none"), p2 + theme(legend.position="none"), nrow=1)
p_legend <- plot_grid(p1_legend, p2_legend, nrow=1)

pdf("Extended_Figure3/image/Extended_Fig3e.pdf", width=7, height=5)
print(p)
dev.off()

pdf("Extended_Figure3/image/Extended_Fig3e_legend.pdf", width=7, height=5)
print(p_legend)
dev.off()








# Violin Plot
options(java.parameters="-Xmx256g")
library(Seurat);library(readr);library(ggplot2);library(cowplot)

Org_1 <- read_rds("Seurat/RData/Org-1.Seu_mm.Epi.SCTransform.SNN.Rds")
Org_2 <- read_rds("Seurat/RData/Org-2.Seu_mm.Epi.SCTransform.SNN.Rds")

Org_1 <- FindClusters(Org_1, resolution = 0.2, verbose=FALSE)
Org_2 <- FindClusters(Org_2, resolution = 0.2, verbose=FALSE)

type_1 <- c("IH", "L", "B", "IL")
Idents(Org_1) <- factor(type_1[Idents(Org_1)], levels=c("B", "IL", "IH", "L"))
type_2 <- c("I", "B", "L", "BP2", "BP1")
Idents(Org_2) <- factor(type_2[Idents(Org_2)], levels=c("BP1", "BP2", "B", "I", "L"))

genes <- c("Krt5", "Krt14", "Trp63", "Krt8", "Krt18", "Upk1a", "Upk1b", "Upk2", "Upk3a", "Krt20")

color_1 <- c("#00BA38", "#61D0FF", "#6166FF", "#F8766D")
color_2 <- c("grey 50", "grey80", "#00BA38", "#619CFF", "#F8766D")

for(i in seq(length(genes))){
  p2 <- VlnPlot(Org_1, features=genes[i], pt.size=0) + scale_fill_manual(values=color_1) + 
          theme(legend.position="none", axis.text.x=element_text(angle=0, hjust=0.5)) + labs(title=NULL, x="7d assembloid")
  p1 <- VlnPlot(Org_2, features=genes[i], pt.size=0) + scale_fill_manual(values=color_2) + 
          theme(legend.position="none", axis.text.x=element_text(angle=0, hjust=0.5)) + labs(title=NULL, x="1d assembloid")
		
  title <- ggdraw() + draw_label(genes[i], fontface="bold", size=17, hjust=0.5)

  assign(genes[i], plot_grid(title, plot_grid(p1, p2, nrow=1), ncol=1, rel_heights=c(0.1,1)))
}

pdf("Figures/Extended Figure 5c violin.pdf", width=24, height=10)
print(plot_grid(Krt5, Krt14, Trp63, Krt8, Krt18, Upk1a, Upk1b, Upk2, Upk3a, Krt20, nrow=3))
dev.off()




































