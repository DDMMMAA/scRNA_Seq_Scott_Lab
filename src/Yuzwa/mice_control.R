library(dplyr)
library(Seurat)
library(patchwork)
library(clusterProfiler)
library(clustree)

#####################################
# include Helper function
source("src/helper_function.R")

# include DEGs/marker utilized in Wang_et_al
source("src/Wang_et_al/marker.R")
#####################################
# raw data directory
target_sample_directory <- "data/raw/Yuzwa/Control-10K_RSEC_MolsPerCell_MEX 1/"

#####################################
# Load the data set
mice_control.data <- Read10X(data.dir = target_sample_directory)
mice_control <- CreateSeuratObject(counts = mice_control.data, 
                                   project = "mice_control", 
                                   min.cells = 3, 
                                   min.features = 200)

####################################
# Quality control(by mitochondrial gene proportion)
mice_control[["percent.mt"]] <- PercentageFeatureSet(mice_control, 
                                                     pattern = "^mt-")

# Visualize QC metrics as a violin plot
VlnPlot(mice_control, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)

# Visualize correlation of nCount_RNA w/ respect to percent.mt/nFeature_RNA
plot1 <- FeatureScatter(mice_control, feature1 = "nCount_RNA", 
                        feature2 = "percent.mt")
plot2 <- FeatureScatter(mice_control, feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA")
plot1 + plot2

# Subset to filter low quality cell
mice_control <- subset(mice_control, subset = nFeature_RNA > 200 & 
                         nFeature_RNA < 2500 & percent.mt < 20)

####################################
# Normalization (by log normalize)
mice_control <- NormalizeData(mice_control, 
                              normalization.method = "LogNormalize", 
                              scale.factor = 10000)

####################################
# Calculate and store top variable gene
mice_control <- FindVariableFeatures(mice_control, selection.method = "vst", 
                                     nfeatures = 2000)

# Identify the top 10 variable gene
top10 <- head(VariableFeatures(mice_control), 10)

# Plot variable features
plot1 <- VariableFeaturePlot(mice_control)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

####################################
#Scaling the data
all.genes <- rownames(mice_control)
mice_control <- ScaleData(mice_control, features = all.genes)

####################################
#(Barely understand the math and graph of this part)
#Linear dimensional reduction
mice_control <- RunPCA(mice_control, 
                       features = VariableFeatures(object = mice_control))

# Visualization of PCA result
VizDimLoadings(mice_control, dims = 1:2, reduction = "pca")
DimPlot(mice_control, reduction = "pca", label = TRUE)
DimHeatmap(mice_control, dims = 1, cells = 500, balanced = TRUE)

###################################
#Determine the 'dimensionality' of the data set
ElbowPlot(mice_control)

###################################
#Cluster the cells
mice_control <- FindNeighbors(mice_control, dims = 1:10)

# Determine optimal resolution
mice_control <- FindClusters.range(mice_control)
clustree(mice_control)

# Set clustering result under resolution = 0.5
mice_control <- FindClusters(mice_control, resolution = 0.5)
head(Idents(mice_control), 5)

###################################
#UMAP!!
mice_control <- RunUMAP(mice_control, dims = 1:10)
DimPlot(mice_control, reduction = "umap", label = TRUE)

##################################
# Finding differentially expressed features(Cluster biomarkers)
mice_control.markers <- FindAllMarkers(mice_control, only.pos = TRUE)

# filter all marker with avg_log2FC > 1
mice_control.markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1)

# Visualization
VlnPlot(mice_control, features = c("Enpp2", "Igfbp2"))

# separate markers by cluster
mice_control.markers.cluster <- split(mice_control.markers, mice_control.markers$cluster)

# subsetting
filtered_deg_list <- lapply(mice_control.markers.cluster, function(df) {
  subset(df, 
         pct.1 > 0.1 & 
           pct.2 > 0.1 & 
           p_val < 0.05 & 
           (avg_log2FC > 0.25 | avg_log2FC < -0.25))
})

# calculate p_FC value
filtered_deg_list <- lapply(filtered_deg_list, function(df) {
  df$p_FC <- (1 - df$p_val) * df$avg_log2FC
  return(df)
})

###################################
# Test of known marker gene
# Astrocyte marker
FeaturePlot(mice_control, features = c("S100b", "Gfap", "Aldh1l1", 
                                       "Sox9", "Gja1", "Hepacam"))

# Oligodendrocyte marker
FeaturePlot(mice_control, features = c("Sox10", "Olig1", "Olig2", "Mbp"))

# Neuron stem marker
FeaturePlot(mice_control, features = c("Sox2", "Vim", "Nes"))

# Neuron marker
FeaturePlot(mice_control, features = c("Neurod2", "Neurod1", "Dcx", "Sp8", 
                                       "Sp9", "Grm8", "Gabra1"))

# Oligo progenitor marker
FeaturePlot(mice_control, features = c("Cspg4"))

# Micro glia marker
FeaturePlot(mice_control, features = c("Ptprc", "Cx3cr1", "Tmem9", "Aif1"))

# Epidermal marker
FeaturePlot(mice_control, features = c("Foxj1"))

##################################
# Gene ontology (Among Top 30 gene of each cluster)
# Run GO for each cluster
GOresult <- enrichGO.clusters(mice_control.markers, 30, "BP")

#################################
# Assign cell identity
new.cluster.ids <- c("Immune_1*", "1", "2", "Oligogendrocyte", "Neuron", "Neuro Stem/OPC",
                     "Cell cycle*", "Astrocyte", "Immune_2*", "Myeloid", "10", "11")
names(new.cluster.ids) <- levels(mice_control)
mice_control <- RenameIdents(mice_control, new.cluster.ids)
DimPlot(mice_control, reduction = "umap", label = TRUE, pt.size = 0.5)
################################
# Test DEGs identified from Wang_et_al Yuzwa data

# Microglia
FeaturePlot(mice_control, features = Wang_et_al_Marker$sup_table3b$Microglia) + 
  plot_annotation(title = "Microglia")

# Astrocyte
FeaturePlot(mice_control, features = Wang_et_al_Marker$sup_table3b$Astro_df) + 
  plot_annotation(title = "Astrocyte_common")

# OPC
FeaturePlot(mice_control, features = Wang_et_al_Marker$sup_table3b$OPC) + 
  plot_annotation(title = "OPC")

# Oligo
FeaturePlot(mice_control, features = Wang_et_al_Marker$sup_table3b$Oligodendrocyte) + 
  plot_annotation(title = "Oligodendrocyte")

# Pre-Oligo
FeaturePlot(mice_control, features = Wang_et_al_Marker$sup_table3b$`Oligodendrocyte-Immature`) + 
  plot_annotation(title = "Pre-oligo")

# Vascular
FeaturePlot(mice_control, features = Wang_et_al_Marker$sup_table3b$Vascular) + 
  plot_annotation(title = "Vascular")

# RG
FeaturePlot(mice_control, features = Wang_et_al_Marker$sup_table3b$RG_df) + 
  plot_annotation(title = "RG-Common")

# Unknown
FeaturePlot(mice_control, features = Wang_et_al_Marker$sup_table3b$Unknown) + 
  plot_annotation(title = "Unknown")

# Tri_IPC
FeaturePlot(mice_control, features = Wang_et_al_Marker$sup_table3b$`IPC-Glia`) + 
  plot_annotation(title = "Tri_IPC")

# CR
FeaturePlot(mice_control, features = Wang_et_al_Marker$sup_table3b$`Cajal-Retzius cell`) + 
  plot_annotation(title = "Cajal-Retzius cell")

# EN
FeaturePlot(mice_control, features = Wang_et_al_Marker$sup_table3b$EN_df) + 
  plot_annotation(title = "EN")

# IN
FeaturePlot(mice_control, features = Wang_et_al_Marker$sup_table3b$IN_df) + 
  plot_annotation(title = "IN")
################################
# Save result seurat object
saveRDS(mice_control, file = "data/processed/Yuzwa/mice_control/seuratObj_mice_control.RData")
