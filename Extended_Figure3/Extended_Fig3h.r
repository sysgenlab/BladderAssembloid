### Extended Fig3 h ###
options(java.parameters="-Xmx256g")
library(Seurat)
library(readr)
library(ggplot2)
library(cowplot)

Seu.Org_day1 <- read_rds("data/scRNA/Seurat/RData/Org_day1.Seu_mm.Epi.SCTransform.SNN.Rds")
Seu.Org_day7 <- read_rds("data/scRNA/Seurat/RData/Org_day7.Seu_mm.Epi.SCTransform.SNN.Rds")

# Rename Clusters
levels(Idents(Seu.Org_day1)) <- c("Basal\nProliferative 1", "Basal\nProliferative 2", "Basal", "Intermediate", "Luminal")
levels(Idents(Seu.Org_day7)) <- c("Basal", "Intermediate\nLow", "Intermediate\nHigh", "Luminal")

# Population Statistics
df <- summary(Idents(Seu.Org_day1))
df1 <- data.frame(x=names(df), y=df/sum(df)*100)
df1$x <- factor(df1$x, levels=df1$x)

df <- summary(Idents(Seu.Org_day7))
df2 <- data.frame(x=names(df), y=df/sum(df)*100)
df2$x <- factor(df2$x, levels=df2$x)

p1 <- ggplot(df1, aes(x=x, y=y)) + geom_bar(stat="identity", width=0.5, fill="royalblue")
p1 <- p1 + theme_classic()
p1 <- p1 + labs(title="1-day bladder assembloid", x=NULL, y="Proportion (%)")
p1 <- p1 + scale_y_continuous(expand=c(0,0), limits=c(0,40))
p1 <- p1 + theme(plot.title=element_text(size=20, hjust=0.5),
                 axis.title=element_text(size=18),
			     axis.text=element_text(size=15))

p2 <- ggplot(df2, aes(x=x, y=y)) + geom_bar(stat="identity", width=0.5, fill="gold")
p2 <- p2 + theme_classic()
p2 <- p2 + labs(title="7-days bladder assembloid", x=NULL, y="Proportion (%)")
p2 <- p2 + scale_y_continuous(expand=c(0,0), limits=c(0,35))
p2 <- p2 + theme(plot.title=element_text(size=20, hjust=0.5),
                 axis.title=element_text(size=18),
			     axis.text=element_text(size=15))

p <- plot_grid(p1, p2, ncol=1)

pdf("Extended_Figure3/image/Extended_Fig3h.pdf", width=7, height=6)
print(p)
dev.off()
