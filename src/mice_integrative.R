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
library(data.table)
library(clustifyr)
library(clustifyrdatahub)

#######################
# include Helper function
source("src/helper_function.R")

#######################
# Read control & treatment Seurat obj result from previous separate analysis
mice_control <- readRDS(file = "data/processed/mice_control/seuratObj_mice_control.RData")
mice_treatment <- readRDS(file = "data/processed/mice_treatment/seuratObj_mice_treatment.RData")

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
ElbowPlot(mice_merged)

mice_merged <- FindNeighbors(mice_merged, 
                             dims = 1:15, 
                             reduction = "pca")
mice_merged <- FindClusters(mice_merged, 
                            resolution = 0.5, 
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

# Determine the optimal resolution
mice_merged <- FindClusters.range(mice_merged, start = 0, end = 1.5, margin = 0.1)
clustree(mice_merged)

# Set clustering result under resolution = 0.5
mice_merged <- FindClusters(mice_merged, resolution = 0.5)
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
            reduction = "umap")

# Oligodendrocyte marker
FeaturePlot(mice_merged, 
            features = c("Sox10", "Olig1", "Olig2", "Mbp"), 
            reduction = "umap")

# Neuron stem marker
FeaturePlot(mice_merged, 
            features = c("Sox2", "Vim", "Nes"), 
            reduction = "umap")

# Neuron marker
FeaturePlot(mice_merged, 
            features = c("Neurod2", "Neurod1", "Dcx", 
                         "Sp8", "Sp9", "Grm8", "Gabra1"), 
            reduction = "umap")

# Oligo progenitor marker
FeaturePlot(mice_merged, 
            features = c("Cspg4", "Pdgfra"), 
            reduction = "umap")

# Micro glia marker
FeaturePlot(mice_merged, 
            features = c("Ptprc", "Cx3cr1", "Tmem9", "Aif1"), 
            reduction = "umap")

# Epidermal marker
FeaturePlot(mice_merged, 
            features = c("Foxj1"), 
            reduction = "umap")

# GO analysis
GOresult <- enrichGO.clusters(mice_merged.markers, 30, "BP")

new.cluster.ids <- c("Proliferating_Microglia*", "1", "Mitochondrial", "Neuron", 
                     "Neuro Stem", "Cell_cycle_1*", "Myeloid", "Cell_cycle_2*", 
                     "Oligo/OPC", "Asterocyte", "Mature_Microglia*", "11", "T-cell*", 
                     "13", "Blood-brain_barrier*", "Cell_movement*")
names(new.cluster.ids) <- levels(mice_merged)
mice_merged <- RenameIdents(mice_merged, new.cluster.ids)
DimPlot(mice_merged, reduction = "umap", label = TRUE, pt.size = 0.5)

# Automated Cell annotation by clustifyr
# using clustering result at res = 0.5 and Mouse Cell Atlas as reference
cor_to_call(clustify(
  input = mice_merged,
  metadata = mice_merged@meta.data$RNA_snn_res.0.5,
  ref_mat = ref_MCA()
))

# Rename cluster based on result above
new.cluster.ids <- c("Macrophage_Klf2_high_1", "Macrophage_Klf2_high_2", 
                     "Macrophage_Klf2_high_3", "Schwann_cell_1", 
                     "Granule_neurons", "Neuron_Kpna2_high", "Schwann_cell_2", 
                     "Proliferating_thymocyte", 
                     "Myelinating_oligodendrocyte", 
                     "Astroglial_cell(Bergman glia)", 
                     "Macrophage_Klf2_high_4", "Astrocyte_Mfe_high", 
                     "Monocyte(Spleen)", 
                     "Schwann_cell_3", "Ductal_cell(Pancreas)", 
                     "Hypothalamic_ependymal cell")
names(new.cluster.ids) <- levels(mice_merged)
mice_merged <- RenameIdents(mice_merged, new.cluster.ids)
DimPlot(mice_merged, reduction = "umap", label = TRUE, pt.size = 0.5)
################################
# Export mice_merged seurat obj into CSV file
data_to_write_out <- as.data.frame(as.matrix(mice_merged[["RNA"]]$data))
fwrite(x = data_to_write_out, row.names = TRUE, file = "mice_merged.csv")

################################
#Identify conserved cell type marker across control & treatment group
t_cell.markers <- FindConservedMarkers(mice_merged, ident.1 = "T-cell*", grouping.var = "orig.ident")
