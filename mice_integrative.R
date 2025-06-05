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
library(clustree)

#######################
# Read control & treatment Seurat obj result from previous separate analysis
mice_control <- readRDS(file = "data/processed/seurat_object/seurat_mice_control.RData")
mice_treatment <- readRDS(file = "data/processed/seurat_object/seurat_mice_treatment.RData")

# Merge two Seurat object
mice_merged <- merge(x = mice_control, 
                     y = mice_treatment)

########################
# Perform analysis w/o integration
# run standard anlaysis workflow (QC is skipped since it's already performed prior merging)
mice_merged <- NormalizeData(mice_merged)
mice_merged <- FindVariableFeatures(mice_merged)
mice_merged <- ScaleData(mice_merged)
mice_merged <- RunPCA(mice_merged)

# Determine the dimantionality
ElbowPlot(mice_control)

mice_merged <- FindNeighbors(mice_merged, 
                             dims = 1:15, 
                             reduction = "pca")
mice_merged <- FindClusters(mice_merged, 
                            resolution = 0.95, 
                            cluster.name = "unintegrated_clusters")
mice_merged <- RunUMAP(mice_merged, 
                       dims = 1:15, 
                       reduction = "pca", 
                       reduction.name = "umap.unintegrated")
DimPlot(mice_merged, 
        reduction = "umap.unintegrated", 
        group.by = c("seurat_clusters"), 
        label = TRUE)
########################
# Perform integration
mice_merged <- IntegrateLayers(object = mice_merged, 
                               method = CCAIntegration, 
                               orig.reduction = "pca", 
                               new.reduction = "integrated.cca", 
                               verbose = TRUE)

# re-join layers after integration
mice_merged[["RNA"]] <- JoinLayers(mice_merged[["RNA"]])

mice_merged <- FindNeighbors(mice_merged, 
                             reduction = "integrated.cca", 
                             dims = 1:15)

# Clustering under different resolution and use clustree package
# to determine an optimal resolution

mice_merged <- FindClusters(mice_merged, 
                            resolution = 1.05)
mice_merged <- RunUMAP(mice_merged, dims = 1:15, reduction = "integrated.cca")
DimPlot(mice_merged, 
        reduction = "umap", 
        split.by = "orig.ident", 
        label = TRUE)

###############################
# Annotation
mice_merged.markers <- FindAllMarkers(mice_merged, only.pos = TRUE)
mice_merged.markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1)

# test known marker gene
# Astrocyte marker
FeaturePlot(mice_merged, 
            features = c("S100b", "Gfap", "Aldh1l1", "Sox9", "Gja1", "Hepacam"), 
            reduction = "umap", 
            split.by = "orig.ident")

# Oligodendrocyte marker
FeaturePlot(mice_merged, 
            features = c("Sox10", "Olig1", "Olig2", "Mbp"), 
            reduction = "umap", 
            split.by = "orig.ident")

# Neuron stem marker
FeaturePlot(mice_merged, 
            features = c("Sox2", "Vim", "Nes"), 
            reduction = "umap", 
            split.by = "orig.ident")

# Neuron marker
FeaturePlot(mice_merged, 
            features = c("Neurod2", "Neurod1", "Dcx", 
                         "Sp8", "Sp9", "Grm8", "Gabra1"), 
            reduction = "umap", 
            split.by = "orig.ident")

# Oligo progenitor marker
FeaturePlot(mice_merged, 
            features = c("Cspg4", "Pdgfra"), 
            reduction = "umap", 
            split.by = "orig.ident")

# Micro glia marker
FeaturePlot(mice_merged, 
            features = c("Ptprc", "Cx3cr1", "Tmem9", "Aif1"), 
            reduction = "umap", 
            split.by = "orig.ident")

# Epidermal marker
FeaturePlot(mice_merged, 
            features = c("Foxj1"), 
            reduction = "umap", 
            split.by = "orig.ident")

# GO analysis
for (i in 0:23) {
  cluster_name <- paste("Cluster", i, sep = "")
  result_name <- paste("result_cluster", i, sep = "")
  assign(cluster_name, head(subset(mice_merged.markers, cluster == i)[["gene"]], 30))
  buf <- get(cluster_name)
  assign(result_name, enrichGO(gene = buf, keyType = "SYMBOL", 
                               OrgDb = "org.Mm.eg.db", ont = "BP"))
  print(i)
  rm(buf)
  rm(cluster_name)
  rm(result_name)
}
rm(i)

#new.cluster.ids <- c("Immune_1*", "1", "Mitochondrial", "Neuron", 
#                     "Neurogenesis*", "Cell Cycle*", "Myeloid", "Cell cycle*", 
#                     "Oligo/OPC", "Asterocyte", "Immune_2*", "11", "T-cell*", 
#                     "13", "Blood-brain barrier*", "Cell movement*")
#names(new.cluster.ids) <- levels(mice_merged)
#mice_merged <- RenameIdents(mice_merged, new.cluster.ids)
DimPlot(mice_merged, reduction = "umap", label = TRUE, pt.size = 0.5)

################################
#Identify conserved cell type marker across control & treatment group
t_cell.markers <- FindConservedMarkers(mice_merged, ident.1 = "T-cell*", grouping.var = "orig.ident")
