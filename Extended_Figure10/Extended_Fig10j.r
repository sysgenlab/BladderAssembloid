### Extended Fig10 j ###
options(stringsAsFactors=FALSE) 
options(java.parameters="-Xmx16g")
library(pheatmap)
library(gplots)

myTPM <- read.delim("data/RNA/TPM.RSEM.DESeq2_normalized.gene.txt", row.names=1)
colnames(myTPM) <- c("L1/FOXA1 k.o tumour", 
                      "B1", 
                      "L1/BMP k.o CAF", 
                      "B1/BMP o.e CAF-FOXA1 k.o tumour", 
                      "L2/FOXA1 k.o tumour", 
                      "B2", 
                      "L2/BMP k.o CAF", 
                      "B2/BMP o.e CAF-FOXA1 k.o tumour", 
                      "L1", 
                      "B1/FOXA1 o.e tumour", 
                      "L1/BMP k.o CAF-FOXA1 o.e tumour", 
                      "B1/BMP o.e CAF", 
                      "L2", 
                      "B2/FOXA1 o.e tumour", 
                      "L2/BMP k.o CAF-FOXA1 o.e tumour", 
                      "B2/BMP o.e CAF")

anno <- read.delim("data/RNA/Anno.RSEM.gene.txt", row.names=1)
LIST <- read.delim("data/RNA/list_MD.BASE47.txt")
rownames(LIST) <- LIST$GENE

#contr <- read.delim(paste0("result/DESeq2/List.DEA.", contrast[i], ".txt"), row.names=1)

# Selecting gene and sample list
anno_LIST <- anno[anno$external_gene_name %in% LIST$GENE,]
Sample <- c("L1/BMP k.o CAF", "L2/BMP k.o CAF", "L1/BMP k.o CAF-FOXA1 o.e tumour", "L2/BMP k.o CAF-FOXA1 o.e tumour")

# Expression Matrix
LIST_TPM <- myTPM[rownames(anno_LIST), Sample]
rownames(LIST_TPM) <- anno_LIST$external_gene_name
LIST_TPM[is.na(LIST_TPM)] <- 0

# Scaling and Normalization
expr <- scale(t(log10(LIST_TPM+1)))
expr[is.nan(expr)] <- 0
MARKER <- LIST[colnames(expr),2]
#DB <- LIST[colnames(expr),3]
#Adjusted_Pvalue <- contr[rownames(anno_LIST),"Adjusted_Pvalue"]
#Adjusted_Pvalue[is.na(Adjusted_Pvalue)] <- 1

#color.v <- c("darkred", "grey")
#breaks.v <- c(-Inf, log10(0.05), 0.1)
#SIG <- character(0)
#for(p in Adjusted_Pvalue){
#  SIG <- c(SIG, color.v[sum(log10(p) >= breaks.v)])
#}

# Setting Annotation and Color
myGene_anno <- data.frame(MARKER, row.names=colnames(expr))
mycolor <- list(MARKER=c(Basal="gold", Luminal="purple", p53like="black", Claudin="grey"))
#mycolor <- list(SIG=c(darkred="darkred", red="red", orange="orange", pink="pink", grey="grey"),
#  MARKER=c(Basal="gold", Luminal="purple", p53like="black", Claudin="grey"),
#  DB=c(BASE47="blue", MD="red"))

pdf("Extended_Figure10/image/Extended_Fig10j.pdf", width=19, height=5)
pheatmap(expr, border_color=FALSE, clustering_distance_cols="euclidean", clustering_method="single", color=bluered(100),
         annotation_col=myGene_anno, annotation_names_col=FALSE, annotation_colors=mycolor, show_colnames=FALSE, fontsize_row=15)
dev.off()
