### Extended Fig10 n ###
options(stringsAsFactors=FALSE)
options(java.parameters="-Xmx16g")
library(ggplot2)
library(ggrepel)

myData <- read.table("data/AME/ame_result_decreased.tsv")

adj_p <- myData$V7[2:nrow(myData)]
score <- numeric(0)
for(j in 1:length(adj_p)){
  tmp <- as.numeric(strsplit(adj_p[j], "e-")[[1]])
  score_tmp <- tmp[2] - log10(tmp[1])
  score <- c(score, score_tmp)
}

colnames(myData) <- myData[1,]
myData <- myData[-1,]
res <- cbind(myData[,c(1,3,4,9,14:17)], score)
rownames(res) <- res$rank
res$rank <- as.numeric(res$rank)

p <- ggplot(res, aes(x=rank, y=score)) + geom_point(size=2)
p <- p + theme_bw()
title <- "Regions with increased\nATAC/H3K27ac peaks in basal"
p <- p + labs(x="JASPAR 2018 (n=1,404)", y=expression("-" ~ log[10] ~ "q-value"), title=title)
p <- p + theme(plot.title=element_text(size=20, hjust=0.5),
               axis.title=element_text(size=18),
			   axis.text.y=element_text(size=15),
               axis.text.x=element_blank(),
			   axis.ticks.x=element_blank())
p <- p + geom_text_repel(data=res[1:10,], 
                         aes(label=res$motif_alt_ID[1:10]), 
                         point.padding=unit(0.3, "lines"), 
                         box.padding=unit(0.8, "lines"))

pdf("Extended_Figure10/image/Extended_Fig10n.pdf", width=5, height=5)
print(p)
dev.off()
