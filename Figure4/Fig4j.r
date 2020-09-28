### Fig 4j ###
options(stringsAsFactors=FALSE)
library(data.table)
library(ggplot2)
library(cowplot)

contr <- "B_BMP_o.e_CAF_vs_B_BMP_o.e_CAF_FOXA1_k.o"
Basal <- "Basal/\nBMP o.e CAF-\nFOXA1 k.o tumour"
Luminal <- "Basal/\nBMP o.e CAF"

ex <- c("ATAC", "H3K27ac", "FOXA1")
title <- c("ATAC", "H3K27ac ChIP", "FOXA1 ChIP")
names(title) <- ex

for(i in seq(ex)){
  # Load matrix
  cmd <- paste0("zcat data/deepTools/matrix.", ex[i], ".", contr, ".GAIN.mat.gz")
  myData <- fread(cmd=cmd, skip=1, header=FALSE, data.table=FALSE)
  
  # Kolmogorov-Smirnov test
  ##x: Basal, y: Luminal
  l <- (ncol(myData)-6)/2
  p.value <- ks.test(x=colMeans(myData[,1:l + 6]), y=colMeans(myData[,1:l + 6 +l]), alternative="greater")$p.value

  df <- data.frame(x=seq(-2000,2000,length.out=l), y=colMeans(myData[,1:(2*l) + 6]), Type=rep(c(Basal, Luminal), each=l))

  p <- ggplot(df, aes(x=x, y=y, colour=Type)) + geom_line(size=2)
  p <- p + theme_bw() 
  p <- p + labs(x=NULL, y=NULL, colour=NULL) 
  p <- p + ggtitle(title[ex[i]], subtitle=paste0("KS-test p-value= ", format(p.value, scientific=TRUE, digits=3)))
  p <- p + scale_x_continuous(expand=c(0,0), breaks=c(-2000,0,2000), labels=c("-2kb", "p.c","2kb"))
  p <- p + scale_y_continuous(expand=c(0,0),limits=c(min(df$y)*0.7,max(df$y)*1.02))
  p <- p + scale_colour_manual(values=c("blue", "red"))
  p <- p + theme(panel.grid=element_blank(), 
                 axis.text=element_text(size=15),  
   	             axis.ticks.length=unit(0.2,"cm"),
	             legend.position="right",
	             legend.text=element_text(size=14),
	             plot.title=element_text(size=22, hjust=0.5),
                 plot.subtitle=element_text(size=15, hjust=0.5),
	             plot.margin=margin(1,1,0,0.5,"cm"))
  
  assign(paste0("p", i), p)
}

legend <- get_legend(p1)
p <- plot_grid(p1 + theme(legend.position="none"), 
               p2 + theme(legend.position="none"), 
               p3 + theme(legend.position="none"), 
               get_legend(p1),
               nrow=1,
               rel_widths=c(1,1,1,0.5))

pdf("Figure4/image/Fig4j.pdf", width=15, height=6)
print(p)
dev.off()
