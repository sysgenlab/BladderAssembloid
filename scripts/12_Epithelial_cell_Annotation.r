options(java.parameters="-Xmx256g")
library(Seurat)
library(readr)

# Org_day1 Epithelial cells
Seu <- read_rds("data/scRNA/Seurat/RData/Org_day1.Seu_mm.Epi.SCTransform.SNN.Rds")
Type <- c("Intermediate (I)",
          "Basal (B)",
          "Luminal (L)",
          "Basal proliferative 2 (BP2)",
          "Basal proliferative 1 (BP1)")
Type.lvl <- c("Basal proliferative 1 (BP1)", 
              "Basal proliferative 2 (BP2)", 
			  "Basal (B)", 
			  "Intermediate (I)", 
			  "Luminal (L)")
				 
Seu@meta.data$Type <- factor(Type[Seu@meta.data$SCT_snn_res.0.2], levels=Type.lvl)
Idents(Seu) <- Seu@meta.data$Type

write_rds(Seu, "data/scRNA/Seurat/RData/Org_day1.Seu_mm.Epi.SCTransform.SNN.Rds")


# Org_day7 Epithelial cells
Seu <- read_rds("data/scRNA/Seurat/RData/Org_day7.Seu_mm.Epi.SCTransform.SNN.Rds")
Type <- c("Intermediate high (IH)",
          "Luminal (L)",
          "Basal (B)",
          "Intermediate low (IL)")
Type.lvl <- c("Basal (B)", 
              "Intermediate low (IL)", 
			  "Intermediate high (IH)", 
			  "Luminal (L)")
				 
Seu@meta.data$Type <- factor(Type[Seu@meta.data$SCT_snn_res.0.2], levels=Type.lvl)
Idents(Seu) <- Seu@meta.data$Type

write_rds(Seu, "data/scRNA/Seurat/RData/Org_day7.Seu_mm.Epi.SCTransform.SNN.Rds")