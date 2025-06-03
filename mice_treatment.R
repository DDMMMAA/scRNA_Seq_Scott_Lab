library(dplyr)
library(Seurat)
library(patchwork)
library(clusterProfiler)

#####################################
#Global Variable
target_sample_directory <- "data/raw/mice/Injury-10K_RSEC_MolsPerCell_MEX 1/"

#####################################
# Load the PBMC data set
mice_treatment.data <- Read10X(data.dir = target_sample_directory)
mice_treatment <- CreateSeuratObject(counts = mice_treatment.data, 
                                      project = "mice_treatment", 
                                      min.cells = 3, 
                                      min.features = 200)

####################################
# Quality control(by mitochondrial gene proportion)
mice_treatment[["percent.mt"]] <- PercentageFeatureSet(mice_treatment, pattern = "^mt-")

# Visualize QC metrics as a violin plot
VlnPlot(mice_treatment, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)

# Visualize correlation of nCount_RNA w/ respect to percent.mt/nFeature_RNA
plot1 <- FeatureScatter(mice_treatment, feature1 = "nCount_RNA", 
                        feature2 = "percent.mt")
plot2 <- FeatureScatter(mice_treatment, feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA")
plot1 + plot2

# Subset to filter low quality cell
mice_treatment <- subset(mice_treatment, subset = nFeature_RNA > 200 & 
                            nFeature_RNA < 2500 & percent.mt < 20)

####################################
# Normalization (by log normalize)
mice_treatment <- NormalizeData(mice_treatment, 
                                 normalization.method = "LogNormalize", 
                                 scale.factor = 10000)

####################################
# Calculate and store top variable gene
mice_treatment <- FindVariableFeatures(mice_treatment, selection.method = "vst", nfeatures = 2000)

# Identify the top 10 variable gene
top10 <- head(VariableFeatures(mice_treatment), 10)

# Plot variable features
plot1 <- VariableFeaturePlot(mice_treatment)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

####################################
#Scaling the data
all.genes <- rownames(mice_treatment)
mice_treatment <- ScaleData(mice_treatment, features = all.genes)

####################################
#(Barely understand the math and graph of this part)
#Linear dimensional reduction
mice_treatment <- RunPCA(mice_treatment, features = VariableFeatures(object = mice_treatment))

# Visualization of PCA result
VizDimLoadings(mice_treatment, dims = 1:2, reduction = "pca")
DimPlot(mice_treatment, reduction = "pca", label = TRUE)
DimHeatmap(mice_treatment, dims = 1, cells = 500, balanced = TRUE)

###################################
#Determine the 'dimensionality' of the data set
ElbowPlot(mice_treatment)

###################################
#Cluster the cells
mice_treatment <- FindNeighbors(mice_treatment, dims = 1:15)
mice_treatment <- FindClusters(mice_treatment, resolution = 0.5)
head(Idents(mice_treatment), 5)

###################################
#UMAP!!
mice_treatment <- RunUMAP(mice_treatment, dims = 1:15)
DimPlot(mice_treatment, reduction = "umap", label = TRUE)

##################################
#Finding differentially expressed features(Cluster biomarkers)
mice_treatment.markers <- FindAllMarkers(mice_treatment, only.pos = TRUE)

# filter all marker with avg_log2FC > 1
mice_treatment.markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1)

# Sta test
cluster0.markers <- FindMarkers(mice_treatment, ident.1 = 0, logfc.threshold = 0.25, 
                                test.use = "roc", only.pos = TRUE)

# Visualization
VlnPlot(mice_treatment, features = c("Enpp2", "Igfbp2"))

###################################
# Test of known marker gene
# Astrocyte marker
FeaturePlot(mice_treatment, features = c("S100b", "Gfap", "Aldh1l1", 
                                       "Sox9", "Gja1", "Hepacam"))

# Oligodendrocyte marker
FeaturePlot(mice_treatment, features = c("Sox10", "Olig1", "Olig2", "Mbp"))

# Neuron stem marker
FeaturePlot(mice_treatment, features = c("Sox2", "Vim", "Nes"))

# Neuron marker
FeaturePlot(mice_treatment, features = c("Neurod2", "Neurod1", "Dcx", "Sp8", 
                                       "Sp9", "Grm8", "Gabra1"))

# Oligo progenitor marker
FeaturePlot(mice_treatment, features = c("Cspg4", "Pdgfra"))

# Micro glia marker
FeaturePlot(mice_treatment, features = c("Ptprc", "Cx3cr1", "Tmem9", "Aif1"))

# Epidermal marker
FeaturePlot(mice_treatment, features = c("Foxj1"))

##################################
# Gene ontology (Among Top 30 gene of each cluster)
# Run GO for each cluster
for (i in 0:14) {
  cluster_name <- paste("Cluster", i, sep = "")
  result_name <- paste("result_cluster", i, sep = "")
  assign(cluster_name, head(subset(mice_treatment.markers, cluster == i)[["gene"]], 30))
  buf <- get(cluster_name)
  assign(result_name, enrichGO(gene = buf, keyType = "SYMBOL", 
                               OrgDb = "org.Mm.eg.db", ont = "BP"))
  print(i)
}

####################################
# Assign cell identity
new.cluster.ids <- c("Immune_1*", "1", "Mitochondrial", "Neuron", "Cell cycle*", "Neuro Stem/OPC",
                     "Myeloid", "7", "8", "Immune_2*", "Astrocyte", "Oligodendrocyte", "T-cell*", "13", "ECM*")
names(new.cluster.ids) <- levels(mice_treatment)
mice_treatment <- RenameIdents(mice_treatment, new.cluster.ids)
DimPlot(mice_treatment, reduction = "umap", label = TRUE, pt.size = 0.5)

################################
# Save result seurat object
saveRDS(mice_treatment, file = "data/processed/seurat_object/seurat_mice_treatment.RData")
