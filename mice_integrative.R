library(Seurat)
library(SeuratData)
library(patchwork)
library(multtest)
library(metap)
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)
library(clusterProfiler)

#######################
# Read control & treatment Seurat obj result from separate analysis
mice_control <- readRDS(file = "data/processed/seurat_object/seurat_mice_control.RData")
mice_treatment <- readRDS(file = "data/processed/seurat_object/seurat_mice_treatment.RData")

# Merge Seurat two object
mice_merged <- merge(x = mice_control, y = mice_treatment)

mice_merged[["RNA"]] <- JoinLayers(mice_merged[["RNA"]])
mice_merged <- FindNeighbors(mice_merged, reduction = "integrated.cca", dims = 1:30)
