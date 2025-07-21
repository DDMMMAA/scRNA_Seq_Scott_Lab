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
library(SingleR)
library(celldex)
library(ggpubr)
library(bit64)
library(enrichplot)
library("org.Mm.eg.db")

#######################
# include Helper function
source("src/helper_function.R")
source("src/Wakui_et_al/Wakui_et_al_helper.R")

#####################################
# Raw data dir
# Loop over Experiment repeat to create correspond Seurat obj
for (exp_repeat in 1:4) {
  # Expression matrix
  target_sample_directory <- paste("data/raw/Wakui_et_al/RNA-sequence/Cortex/Cortex_", exp_repeat, sep = "")
  # Barcode group (Whitelist)
  barcode_dir <- paste("data/raw/Wakui_et_al/Whitelist/Cortex/Cortex_", exp_repeat, ".whitelist", sep = "")
  # Load the data set
  assign(paste("Wakui_cortex", ".data", sep = ""), Read10X(data.dir = target_sample_directory))
  assign(paste("Wakui_cortex", exp_repeat, sep = ""), CreateSeuratObject(counts = Wakui_cortex.data, 
                                                                         project = paste("Wakui_cortex", exp_repeat, sep = "")))
  # Align read's barcode w/ treatment group
  assign(paste("Wakui_cortex", exp_repeat, sep = ""), Align_barcode_group(get(paste("Wakui_cortex", exp_repeat, sep = "")), barcode_dir, exp_repeat))
  message(exp_repeat)
}

#######################
# Merge four Seurat object
Wakui_merged <- merge(x = Wakui_cortex1, 
                     y = c(Wakui_cortex2, Wakui_cortex3, Wakui_cortex4))

# QC
# This dataset is seemingly poor, only ~10% cell remained
Wakui_merged <- subset(Wakui_merged, subset = nFeature_RNA > 200 & 
                         nFeature_RNA < 6000)
########################
# Perform analysis w/o integration
# run standard anlaysis workflow
Wakui_merged <- NormalizeData(Wakui_merged)
Wakui_merged <- FindVariableFeatures(Wakui_merged)
Wakui_merged <- ScaleData(Wakui_merged)
Wakui_merged <- RunPCA(Wakui_merged)

# Determine the dimentionality
ElbowPlot(Wakui_merged)

# Dimentionality = 10
Wakui_merged <- FindNeighbors(Wakui_merged, 
                             dims = 1:10, 
                             reduction = "pca")

Wakui_merged <- FindClusters(Wakui_merged,
                             resolution = 1,
                             cluster.name = "unintegrated_clusters")

Wakui_merged <- RunUMAP(Wakui_merged, 
                       dims = 1:10, 
                       reduction = "pca", 
                       reduction.name = "umap.unintegrated")

DimPlot(Wakui_merged, 
        reduction = "umap.unintegrated", 
        split.by = c("orig.ident"), 
        label = TRUE)
########################
# Perform integration
Wakui_merged <- IntegrateLayers(object = Wakui_merged, 
                               method = CCAIntegration, 
                               orig.reduction = "pca", 
                               new.reduction = "integrated.cca", 
                               verbose = TRUE)

# re-join layers after integration
Wakui_merged[["RNA"]] <- JoinLayers(Wakui_merged[["RNA"]])

Wakui_merged <- FindNeighbors(Wakui_merged, 
                             reduction = "integrated.cca", 
                             dims = 1:10)

# Determine the optimal resolution
Wakui_merged <- FindClusters.range(Wakui_merged, start = 0, end = 1.5, margin = 0.1)
clustree(Wakui_merged)

# Set clustering result under resolution = 1
Wakui_merged <- FindClusters(Wakui_merged, resolution = 1)
Wakui_merged <- RunUMAP(Wakui_merged, dims = 1:10, reduction = "integrated.cca")
DimPlot(Wakui_merged, 
        reduction = "umap", 
        label = TRUE)

###############################
# Annotation
Wakui_merged.markers <- FindAllMarkers(Wakui_merged, only.pos = TRUE)
Wakui_merged.markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1)

# test known marker gene
# Astrocyte marker
FeaturePlot(Wakui_merged, 
            features = c("S100b", "Gfap", "Aldh1l1", "Sox9", "Gja1", "Hepacam"), 
            reduction = "umap")

# Oligodendrocyte marker
FeaturePlot(Wakui_merged, 
            features = c("Sox10", "Olig1", "Olig2", "Mbp"), 
            reduction = "umap")

# Neuron stem marker
FeaturePlot(Wakui_merged, 
            features = c("Sox2", "Vim", "Nes"), 
            reduction = "umap")

# Neuron marker
FeaturePlot(Wakui_merged, 
            features = c("Neurod2", "Neurod1", "Dcx", 
                         "Sp8", "Sp9", "Grm8", "Gabra1"), 
            reduction = "umap")

# Oligo progenitor marker
FeaturePlot(Wakui_merged, 
            features = c("Cspg4", "Pdgfra"), 
            reduction = "umap")

# Micro glia marker
FeaturePlot(Wakui_merged, 
            features = c("Ptprc", "Cx3cr1", "Tmem9", "Aif1"), 
            reduction = "umap")

# Epidermal marker
FeaturePlot(Wakui_merged, 
            features = c("Foxj1"), 
            reduction = "umap")

# GO analysis
GOresult <- enrichGO.clusters(Wakui_merged.markers, 30, "BP")

# Manual Annotation
new.cluster.ids <- c("Microglia_1", "Neuron_1", "Cardiac muscle*", "OPC", 
                     "Neuron_2", "Neuron_3", "Neuron_4", "Neuron_5", 
                     "Neuron_6", "Neuron_7", "Cytoskeleton*", "Neuron_8", 
                     "Neuron_9", "Neuron_10", "Neuron_11", "Neuron_12", 
                     "Oligo", "Astro_1", "ECM", "Neuron_13", "Neuron_14", 
                     "Astro_2", "Pre_oligo_1", "Pre_oligo_2", "Pre_oligo_3", 
                     "Fibroblast", "Neuron_15", "Pre_oligo_4", "Microglia_2", 
                     "Neuron_16", "Neuron_17")
names(new.cluster.ids) <- levels(Wakui_merged)
Wakui_merged <- RenameIdents(Wakui_merged, new.cluster.ids)
DimPlot(Wakui_merged, reduction = "umap", label = TRUE, pt.size = 0.5)

# Automated Cell annotation by clustifyr
# using clustering result at res = 0.5 and Mouse Cell Atlas as reference
x <- cor_to_call(clustify(
  input = Wakui_merged,
  metadata = Wakui_merged@meta.data$RNA_snn_res.1,
  ref_mat = ref_MCA()
))

# Rename cluster based on result x above
new.cluster.ids <- c("Macrophage_1", "Dopaminergic neuron_1", "Granule neurons_1", 
                     "Oligodendrocyte precursor_1", "Dopaminergic neurons_2", 
                     "Dopaminergic neurons_3", "Dopaminergic neurons_4", 
                     "Dopaminergic neurons_5", "Dopaminergic neurons_6", 
                     "Granule neuron_2", "Dopaminergic neurons_7", 
                     "Granule neuron_3", "Granule neurons_4", 
                     "Dopaminergic neurons_8", "Postmitotic neurons", 
                     "Hippocampus neurons", "Myelinating oligodendrocyte", 
                     "Astrocyte_1", "Dopaminergic neurons_9", 
                     "Neural progenitor cell_1", "Granule neurons_5", 
                     "Astrocyte_2", "Oligodendrocyte precursor_2", 
                     "Oligodendrocyte precursor_3", "Oligodendrocyte precursor_4", 
                     "Stromal cell", "Granule neurons_6", 
                     "Oligodendrocyte precursor_5", "Macrophage_2", 
                     "Neural progenitor cell_2", "Granule neurons_7")
names(new.cluster.ids) <- levels(Wakui_merged)
Wakui_merged <- RenameIdents(Wakui_merged, new.cluster.ids)
DimPlot(Wakui_merged, reduction = "umap", label = TRUE, pt.size = 0.5)

# Automated Cell annotation by SingleR
ref_mice <- fetchReference("mouse_rnaseq", "2024-02-26")
SingleR_result <- SingleR(as.data.frame(as.matrix(Wakui_merged[["RNA"]]$data)), 
                          ref_mice, ref_mice$label.main)
Wakui_merged$SingleR.labels <- SingleR_result$labels
DimPlot(Wakui_merged, reduction = "umap", label = TRUE, pt.size = 0.5, 
        group.by = "SingleR.labels", repel = T)
################################
# plot cluster proportion
pt <- table(Wakui_merged$seurat_clusters, Wakui_merged$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion")
################################
# Identify conserved cell type marker across control & treatment group
# This step meant to help manual annotation
Idents(Wakui_merged) <- "seurat_clusters"

Wakui_merged$seurat_clusters <- factor(Wakui_merged$seurat_clusters, 
                                      levels = levels(Wakui_merged))

markers.to.plot <- c("S100a11")
DotPlot(Wakui_merged, features = markers.to.plot, cols = c("blue", "red", "yellow", "black"), 
        dot.scale = 8, split.by = "orig.ident") + RotatedAxis()

################################
# Identify DEG across conditions
DEGs_result <- DEGs_across_condition(Wakui_merged, "Sham")

# Subset DEGs that have either pct.1 or pct.2 > 0.1, meant to accommodate low dataset quality
DEGs_result_subset <- subset_DEGs(DEGs_result, 
                                   "pct.1 > 0.1 & 
                                   pct.2 > 0.1 & 
                                   p_val < 0.05 & 
                                   (avg_log2FC > 0.25 | avg_log2FC < -0.25)")

DEGs_result_subset_Up <- subset_DEGs(DEGs_result_subset, 
                                     "avg_log2FC < 0")

DEGs_result_subset_Down <- subset_DEGs(DEGs_result_subset, 
                                     "avg_log2FC > 0")

# Subset Seurat obj so that it only contain "Sham" and "Hypoxia" group
Wakui_merged_hypoxia <- subset(Wakui_merged, subset = orig.ident == "Hypoxia" | orig.ident == "Sham")
# Subset Seurat obj so that it only contain "Sham" and "Hypoxia" group
Wakui_merged_HIE <- subset(Wakui_merged, subset = orig.ident == "HIE" | orig.ident == "Sham")

# Extract top "num_gene" Upregulated DEGs of cluster "cluster_ident", sorted by given "col_name"
# Or can use two lines of commented code below to itreatively GO DEGs of all cluster among all treatment
# DEGs_GO_Up <- DEGs_enrichGO(DEGs_result_subset_Up, export_figure = F, dir = '')
# DEGs_GO_Down <- DEGs_enrichGO(DEGs_result_subset_Down, export_figure = F, dir = '')
genes.to.label <- extract_top_DEGs(DEGs_result_subset_Up, 
                                   cluster_ident = 0, 
                                   num_gene = 100, 
                                   col_name = "p_FC", 
                                   decreasing = F, 
                                   "Hypoxia")

# GO analysis of top DEGs selected above
DEGs_GO <- enrichGO(gene = genes.to.label, 
                    keyType = "SYMBOL", 
                    OrgDb = "org.Mm.eg.db", 
                    ont = "BP")

barplot(DEGs_GO, showCategory = 5)

# DEG visulization
plots <- VlnPlot(Wakui_merged, features = genes.to.label[0:5], 
                 split.by = "orig.ident", group.by = "seurat_clusters",
                 pt.size = 0.1, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

# extract top 5 up & down regulated DEGs based on p_FC value
Selected_DEGs_Up <- list()
k <- 0
for (i in names(DEGs_result_subset_Up)) {
  Selected_DEGs_Up[[i]] <- list()
  for (j in names(DEGs_result_subset_Up[[i]])) {
    Selected_DEGs_Up[[i]][[j]] <- extract_top_DEGs(DEGs_result_subset_Up, 
                                                   cluster_ident = k, 
                                                   num_gene = 5, 
                                                   col_name = "p_FC", 
                                                   decreasing = F, 
                                                   treatment_group = j)
  }
  k <- k + 1
}

Selected_DEGs_Down <- list()
k <- 0
for (i in names(DEGs_result_subset_Down)) {
  Selected_DEGs_Down[[i]] <- list()
  for (j in names(DEGs_result_subset_Down[[i]])) {
    Selected_DEGs_Down[[i]][[j]] <- extract_top_DEGs(DEGs_result_subset_Down, 
                                                   cluster_ident = k, 
                                                   num_gene = 5, 
                                                   col_name = "p_FC", 
                                                   decreasing = T, 
                                                   treatment_group = j)
  }
  k <- k + 1
}

# Plot spreadsheet of top 5 Upregulated DEGs
Selected_DEGs_Up <- as.data.frame(do.call(rbind, Selected_DEGs_Up))
Selected_DEGs_Up <- data.frame(
  Cluster = sub("^DEG_result_", "", rownames(Selected_DEGs_Up)),
  Hypoxia = sapply(Selected_DEGs_Up$Hypoxia, function(x) paste(x, collapse = ", ")),
  CI = sapply(Selected_DEGs_Up$CI, function(x) paste(x, collapse = ", ")), 
  HIE = sapply(Selected_DEGs_Up$HIE, function(x) paste(x, collapse = ", ")),
  stringsAsFactors = FALSE
)
ggtexttable(Selected_DEGs_Up, rows = NULL, theme = ttheme("light"))

# Plot spreadsheet of top 5 Downregulated DEGs
Selected_DEGs_Down <- as.data.frame(do.call(rbind, Selected_DEGs_Down))
Selected_DEGs_Down <- data.frame(
  Cluster = sub("^DEG_result_", "", rownames(Selected_DEGs_Down)),
  Hypoxia = sapply(Selected_DEGs_Down$Hypoxia, function(x) paste(x, collapse = ", ")),
  CI = sapply(Selected_DEGs_Down$CI, function(x) paste(x, collapse = ", ")), 
  HIE = sapply(Selected_DEGs_Down$HIE, function(x) paste(x, collapse = ", ")),
  stringsAsFactors = FALSE
)
ggtexttable(Selected_DEGs_Down, rows = NULL, theme = ttheme("light"))
################################
# Export Wakui_merged seurat obj into CSV file
data_to_write_out <- as.data.frame(as.matrix(Wakui_merged[["RNA"]]$data))
fwrite(x = data_to_write_out, row.names = TRUE, file = "Wakui_merged.csv")