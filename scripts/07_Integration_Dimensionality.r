# tSNE after SCTransform
## Seurat-3.1.0
options(java.parameters="-Xmx32g")
library(Seurat)
library(scater)
library(readr)
library(cowplot)
library(dplyr)

# Load Data
sce <- read_rds("data/scRNA/scater/RData_filt/WT.sce_mm.QC.norm.filt.Rds")
rownames(sce) <- rowData(sce)$symbol
WT <- as.Seurat(sce)
Anno_WT <- rowData(sce)[,1:10]

sce <- read_rds("data/scRNA/scater/RData_filt/Org_day1.sce_mm.QC.norm.filt.Rds")
rownames(sce) <- rowData(sce)$symbol
Org_day1 <- as.Seurat(sce)
Anno_Org1 <- rowData(sce)[,1:10]

sce <- read_rds("data/scRNA/scater/RData_filt/Org_day7.sce_mm.QC.norm.filt.Rds")
rownames(sce) <- rowData(sce)$symbol
Org_day7 <- as.Seurat(sce)
Anno_Org7 <- rowData(sce)[,1:10]

# Build Annotation
Anno <- merge(Anno_WT, Anno_Org1, by=colnames(Anno_WT), all=TRUE)
Anno <- merge(Anno, Anno_Org7, by=colnames(Anno), all=TRUE)
Anno$symbol[duplicated(Anno$symbol)] <- paste0(Anno$symbol[duplicated(Anno$symbol)], "_dup")
rownames(Anno) <- Anno$symbol

Anno <- Anno[,c(4,5,6,8,9,10)]
  
file_Anno <- "data/scRNA/Anno/Anno.Integrated.biomaRt.txt"
write.table(Anno, file_Anno, sep="\t", quote=F, row.names=TRUE, col.names=TRUE)

colnames(WT@meta.data)[which(colnames(WT@meta.data) == "total_counts")] <- "nCount_RNA"
colnames(WT@meta.data)[which(colnames(WT@meta.data) == "total_features_by_counts")] <- "nFeature_RNA"

colnames(Org_day1@meta.data)[which(colnames(Org_day1@meta.data) == "total_counts")] <- "nCount_RNA"
colnames(Org_day1@meta.data)[which(colnames(Org_day1@meta.data) == "total_features_by_counts")] <- "nFeature_RNA"

colnames(Org_day7@meta.data)[which(colnames(Org_day7@meta.data) == "total_counts")] <- "nCount_RNA"
colnames(Org_day7@meta.data)[which(colnames(Org_day7@meta.data) == "total_features_by_counts")] <- "nFeature_RNA"


# Calculate Mitochonrial gene proportions
WT <- PercentageFeatureSet(WT, pattern = "^mt-", col.name="percent.mt")
Org_day1 <- PercentageFeatureSet(Org_day1, pattern = "^mt-", col.name="percent.mt")
Org_day7 <- PercentageFeatureSet(Org_day7, pattern = "^mt-", col.name="percent.mt")

# SCTransform
WT <- SCTransform(WT, vars.to.regress = "percent.mt")
Org_day1 <- SCTransform(Org_day1, vars.to.regress = "percent.mt")
Org_day7 <- SCTransform(Org_day7, vars.to.regress = "percent.mt")

# Integration
sample.list <- c(WT, Org_day1, Org_day7)

sample.features <- SelectIntegrationFeatures(object.list = sample.list, nfeatures = 3000)
sample.list.pre <- PrepSCTIntegration(object.list = sample.list, anchor.features = sample.features)
sample.anchors <- FindIntegrationAnchors(object.list = sample.list.pre, normalization.method = "SCT", anchor.features = sample.features)
integrated <- IntegrateData(anchorset = sample.anchors, normalization.method = "SCT")

# PCA
integrated <- RunPCA(object = integrated, verbose=FALSE)

# Exploring Optimal Clustering Number
integrated <- JackStraw(integrated, num.replicate=100, verbose=FALSE)
integrated <- ScoreJackStraw(integrated, dims = 1:20)

p <- plot_grid(JackStrawPlot(integrated, dims = 1:20), ElbowPlot(integrated))
pdf("data/scRNA/Seurat/Dimensionality/WT-Org.mm.JackStrawPlot.ElbowPlot.pdf", height=5, width=10)
print(p)
dev.off()

# Saving Data
write_rds(integrated, "data/scRNA/Seurat/RData/WT-Org.integrated.Seu_mm.SCTransform.Rds")
