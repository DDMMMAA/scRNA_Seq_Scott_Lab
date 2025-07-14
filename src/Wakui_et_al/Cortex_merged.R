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
Wakui_merged <- subset(Wakui_merged, subset = nFeature_RNA > 200 & 
                         nFeature_RNA < 6000)
########################
# Perform analysis w/o integration
# run standard anlaysis workflow (QC is skipped since it's already performed prior merging)
Wakui_merged <- NormalizeData(Wakui_merged)
Wakui_merged <- FindVariableFeatures(Wakui_merged)
Wakui_merged <- ScaleData(Wakui_merged)
Wakui_merged <- RunPCA(Wakui_merged)

# Determine the dimentionality
ElbowPlot(Wakui_merged)

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

# Rename cluster based on result above
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
DEGs_result <- DEGs_across_condition(Wakui_merged)

# DEG visulization
# Extract top 5 DEGs of each cluster, sorted by increasing P-val
genes.to.label <- extract_top_DEGs(DEGs_result, 
                                   cluster_ident = 0, 
                                   num_gene = 5, 
                                   col_name = "p_val", 
                                   decreasing = F)
plots <- VlnPlot(Wakui_merged, features = genes.to.label, 
                 split.by = "orig.ident", group.by = "seurat_clusters",
                 pt.size = 0.1, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

# Construct DEGs dataframe of manually selected "interesting" DEGs
# Yeah, this is dump, I hope there is an elegant automatic way to do this
Selected_DEGs <- list(Cluster0 = list(Up = c("Cp"), 
                                       Down = c("Msh5")), 
                          cluster1 = list(Up = c("Cp"), 
                                       Down = c("Plac8", "Msh5")), 
                          cluster2 = list(Up = c("Cp"), 
                                       Down = c("Msh5", "Tnfsf18")), 
                          cluster3 = list(Up = c("Xrccl", "Hbb-bt", "Tmsb4x"), 
                                       Down = c("Gm6260", "1500004A13Rik")), 
                          cluster4 = list(Up = c("Rbfox1", "Tatdn1"), 
                                       Down = c("Cwc22", "Acaa2")), 
                          cluster5 = list(Up = c("Cp"), 
                                       Down = c("Gm42047")), 
                          cluster6 = list(Up = c("Cp", "Rbfox1", "Ube216"), 
                                       Down = c("Msh5")), 
                          cluster7 = list(Up = c("Cp", "Ntm", "Lsamp", "Meg3"), 
                                       Down = c("Msh5", "Hpgd", "Rrp9")), 
                          cluster8 = list(Up = c("Hbb-bt", "Hbb-bs", "Hba-a", "Srsf11", "Atpif1", "Lin7a", "Caskin2", "Dele1"), 
                                       Down = c("Brcc3", "Rpl9-ps6", "Fcor", "Stk32c", "Gm5087", "Gm3445")), 
                          cluster9 = list(Up = c("Col1a1", "Col12a1", "Col3a1", "Slc47a1", "Slc22a8", "Slc4a10", "Mgp", "Airn", "Aldh1a2"), 
                                       Down = c("Hbb-bt", "Hbb-bs", "Hba-a")), 
                          cluster10 = list(Up = c("Slc13a3", "Gm15675", "Fpr2", "Ttn"), 
                                        Down = c("Gm49144", "Slc27a6", "Pdgfc", "Ifi27l2a", "Ccl5")), 
                          cluster11 = list(Up = c("Jund", "Veph1", "Arhgap23", "Rps4x", "Cx3cr1"), 
                                        Down = c("Ppp5c", "Zfp335os", "Zfp687", "Pip4p1", "Atg4b")), 
                          cluster12 = list(Up =c("Chil3", "Hba-a1", "Hbb-bt", "Cxcr6", "Trgc2"), 
                                        Down = c("Vpreb1", "Dntt", "Cwc22", "Igll1", "Cd46")), 
                          cluster13 = list(Up = c("E230029C05Rik", "A330076H08Rik", "Cx3cr1", "Mctp1", "Sorcs1"), 
                                        Down = c("Pdcd4", "Rc3h2", "Galnt11", "Sp3", "Ube2g1")), 
                          cluster14 = list(Up = c("Hba-a2"), 
                                        Down = c("Sdhb", "Hbp1", "Rhbdf2")), 
                          cluster15 = list(Up = c("Cx3cr1", "Pla2g7", "Ccser1", "Gabrb1"), 
                                        Down = c("Sncaip", "Dhcr24", "Zfp950", "Cdip1", "Ndufa13"))
)
Selected_DEGs <- as.data.frame(do.call(rbind, Selected_DEGs))
Selected_DEGs <- data.frame(
  Cluster = rownames(Selected_DEGs),
  Up = sapply(Selected_DEGs$Up, function(x) paste(x, collapse = ", ")),
  Down = sapply(Selected_DEGs$Down, function(x) paste(x, collapse = ", ")),
  stringsAsFactors = FALSE
)
ggtexttable(Selected_DEGs, rows = NULL, theme = ttheme("light"))
################################
# Export Wakui_merged seurat obj into CSV file
data_to_write_out <- as.data.frame(as.matrix(Wakui_merged[["RNA"]]$data))
fwrite(x = data_to_write_out, row.names = TRUE, file = "Wakui_merged.csv")