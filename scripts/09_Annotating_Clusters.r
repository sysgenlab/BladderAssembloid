options(java.parameters="-Xmx256g")
library(Seurat)
library(readr)
library(dplyr)
library(xlsx)

# WT
Seu <- read_rds("data/scRNA/Seurat/RData/WT.Seu_mm.SCTransform.SNN.Rds")
Type <- c("Epithelial cells", 
          "Stromal fibroblasts", 
		  "Smooth muscle cells", 
		  "Endothelial cells", 
		  "Schwann cells", 
		  "Endothelial cells", 
		  "Smooth muscle cells", 
		  "Macrophages", 
		  "Mesothelial cells")
Subtype <- c("Epithelial cells", 
             "Stromal fibroblasts", 
			 "Smooth muscle cells 1", 
			 "Endothelial cells 1", 
			 "Schwann cells", 
			 "Endothelial cells 2", 
			 "Smooth muscle cells 2", 
			 "Macrophages", 
			 "Mesothelial cells")
Type.lvl <- c("Epithelial cells", 
              "Stromal fibroblasts", 
			  "Smooth muscle cells", 
			  "Endothelial cells", 
			  "Schwann cells", 
			  "Macrophages", 
			  "Mesothelial cells")
Subtype.lvl <- c("Epithelial cells", 
                 "Stromal fibroblasts", 
				 "Smooth muscle cells 1", 
				 "Smooth muscle cells 2", 
				 "Endothelial cells 1", 
				 "Endothelial cells 2", 
				 "Schwann cells", 
				 "Macrophages", 
				 "Mesothelial cells")

Seu@meta.data$Type <- factor(Type[Seu@meta.data$SCT_snn_res.0.1], levels=Type.lvl)
Seu@meta.data$Subtype <- factor(Subtype[Seu@meta.data$SCT_snn_res.0.1], levels=Subtype.lvl)
Idents(Seu) <- Seu@meta.data$Type

write_rds(Seu, "data/scRNA/Seurat/RData/WT.Seu_mm.SCTransform.SNN.Rds")

## Find Marker Genes for each Type (Wilcox)
Anno <- read.delim("data/scRNA/Anno/Anno.WT.mm.biomaRt.txt", row.names=1)
markers <- FindAllMarkers(Seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose=FALSE)
MG <- markers %>% group_by(cluster)
res <- cbind(Anno[MG$gene, ], as.data.frame(MG))
res <- res[res$p_val_adj < 0.05,]
res <- split(res, f=res$cluster)
  
file_MG <- "data/scRNA/Seurat/MG_list/Marker.Genes.Wilcox.WT.mm.Type.xlsx"
for(c in seq(length(res))){
  Type <- as.character(unique(res[[c]]$cluster))
  write.xlsx2(res[[c]], file_MG, row.names=FALSE, append=TRUE, sheetName=Type)
}


# Org_day1
Seu <- read_rds("data/scRNA/Seurat/RData/Org_day1.Seu_mm.SCTransform.SNN.Rds")
Type <- c("Smooth muscle cells",
          "Stromal fibroblasts",
		  "Epithelial cells",
		  "Schwann cells",
		  "Endothelial cells",
		  "Others (cycling cells)",
		  "Macrophages",
          "Endothelial cells")
Type.lvl <- c("Epithelial cells", 
              "Stromal fibroblasts", 
			  "Smooth muscle cells", 
			  "Endothelial cells", 
			  "Schwann cells", 
              "Macrophages",
			  "Others (cycling cells)")
Subtype <- c("Smooth muscle cells",
             "Stromal fibroblasts",
			 "Epithelial cells",
			 "Schwann cells",
			 "Endothelial cells",
			 "Others (cycling cells)", 
             "Macrophages",
             "Endothelial cells")
Subtype.lvl <- c("Epithelial cells",
				 "Stromal fibroblasts",
				 "Smooth muscle cells",
				 "Endothelial cells",
				 "Schwann cells",
				 "Macrophages",
				 "Others (cycling cells)")
				 
Seu@meta.data$Type <- factor(Type[Seu@meta.data$SCT_snn_res.0.1], levels=Type.lvl)
Seu@meta.data$Subtype <- factor(Subtype[Seu@meta.data$SCT_snn_res.0.1], levels=Subtype.lvl)
Idents(Seu) <- Seu@meta.data$Type

write_rds(Seu, "data/scRNA/Seurat/RData/Org_day1.Seu_mm.SCTransform.SNN.Rds")

## HULEC
Seu <- read_rds("data/scRNA/Seurat/RData/Org_day1.Seu_hs.SCTransform.SNN.Rds")
Type <- "HULEC"
Seu@meta.data$Type <- factor(Type)
Idents(Seu) <- Seu@meta.data$Type

write_rds(Seu, "data/scRNA/Seurat/RData/Org_day1.Seu_hs.SCTransform.SNN.Rds")


# Org_day7
Seu <- read_rds("data/scRNA/Seurat/RData/Org_day7.Seu_mm.SCTransform.SNN.Rds")
Type <- c("Smooth muscle cells", 
          "Stromal fibroblasts", 
		  "Epithelial cells", 
		  "Epithelial cells", 
		  "Stromal fibroblasts", 
		  "Smooth muscle cells", 
		  "Endothelial cells", 
		  "Endothelial cells", 
		  "Schwann cells", 
		  "Macrophages")
Type.lvl <- c("Epithelial cells", 
              "Stromal fibroblasts", 
			  "Smooth muscle cells", 
			  "Endothelial cells", 
			  "Schwann cells", 
			  "Macrophages")
Subtype <- c("Smooth muscle cells 1", 
             "Stromal fibroblasts 1", 
			 "Epithelial cells 1", 
			 "Epithelial cells 2", 
			 "Stromal fibroblasts 2", 
			 "Smooth muscle cells 2", 
			 "Endothelial cells 1", 
			 "Endothelial cells 2", 
			 "Schwann cells", 
			 "Macrophages")
Subtype.lvl <- c("Epithelial cells 1",
                 "Epithelial cells 2",
				 "Stromal fibroblasts 1",
				 "Stromal fibroblasts 2",
				 "Smooth muscle cells 1",
				 "Smooth muscle cells 2",
				 "Endothelial cells 1",
				 "Endothelial cells 2",
				 "Schwann cells",
				 "Macrophages")
Seu@meta.data$Type <- factor(Type[Seu@meta.data$SCT_snn_res.0.1], levels=Type.lvl)
Seu@meta.data$Subtype <- factor(Subtype[Seu@meta.data$SCT_snn_res.0.1], levels=Subtype.lvl)
Idents(Seu) <- Seu@meta.data$Type

write_rds(Seu, "data/scRNA/Seurat/RData/Org_day7.Seu_mm.SCTransform.SNN.Rds")

## HULEC
Seu <- read_rds("data/scRNA/Seurat/RData/Org_day7.Seu_hs.SCTransform.SNN.Rds")
Type <- "HULEC"
Seu@meta.data$Type <- factor(Type)
Idents(Seu) <- Seu@meta.data$Type

write_rds(Seu, "data/scRNA/Seurat/RData/Org_day7.Seu_hs.SCTransform.SNN.Rds")


# WT-Org Integrated
Seu <- read_rds("data/scRNA/Seurat/RData/WT-Org.integrated.Seu_mm.SCTransform.SNN.Rds")
Type <- c("Stromal fibroblasts",
          "Epithelial cells",
		  "Smooth muscle cells",
		  "Endothelial cells",
		  "Smooth muscle cells",
		  "Schwann cells",
		  "Macrophages")
Type.lvl <- c("Epithelial cells",
              "Stromal fibroblasts",
			  "Smooth muscle cells",
			  "Endothelial cells",
			  "Schwann cells",
			  "Macrophages")
Subtype <- c("Stromal fibroblasts",
             "Epithelial cells",
	  	     "Smooth muscle cells 1",
		     "Endothelial cells",
		     "Smooth muscle cells 2",
		     "Schwann cells",
		     "Macrophages")
Subtype.lvl <- c("Epithelial cells",
                 "Stromal fibroblasts",
			     "Smooth muscle cells 1",
				 "Smooth muscle cells 2",
			     "Endothelial cells",
			     "Schwann cells",
			     "Macrophages")

Seu@meta.data$Type <- factor(Type[Seu@meta.data$integrated_snn_res.0.1], levels=Type.lvl)
Seu@meta.data$Subtype <- factor(Subtype[Seu@meta.data$integrated_snn_res.0.1], levels=Subtype.lvl)
Idents(Seu) <- Seu@meta.data$Type

write_rds(Seu, "data/scRNA/Seurat/RData/WT-Org.integrated.Seu_mm.SCTransform.SNN.Rds")
