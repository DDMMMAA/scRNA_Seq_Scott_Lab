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

#######################
# include Helper function
source("src/helper_function.R")

# include DEGs/marker utilized in Wang_et_al
source("src/Wang_et_al/marker.R")
#######################
# Read control & treatment Seurat obj result from previous separate analysis
mice_control <- readRDS(file = "data/processed/Yuzwa/mice_control/seuratObj_mice_control.RData")
mice_treatment <- readRDS(file = "data/processed/Yuzwa/mice_treatment/seuratObj_mice_treatment.RData")

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

mice_merged.markers.cluster <- split(mice_merged.markers, mice_merged.markers$cluster)

filtered_deg_list <- lapply(mice_merged.markers.cluster, function(df) {
  subset(df, 
         pct.1 > 0.1 & 
           pct.2 > 0.1 & 
           p_val < 0.05 & 
           (avg_log2FC > 0.25 | avg_log2FC < -0.25))
})

filtered_deg_list <- lapply(filtered_deg_list, function(df) {
  df$p_FC <- (1 - df$p_val) * df$avg_log2FC
  return(df)
})

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
FeaturePlot(mice_control, 
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

# Manual Annotation
new.cluster.ids <- c("Proliferating_Microglia*", "1", "Mitochondrial", "Neuron", 
                     "Neuro Stem", "Cell_cycle_1*", "Myeloid", "Cell_cycle_2*", 
                     "Oligo/OPC", "Asterocyte", "Mature_Microglia*", "11", 
                     "T-cell*", 
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

# Automated Cell annotation by SingleR
ref_mice <- fetchReference("mouse_rnaseq", "2024-02-26")
SingleR_result <- SingleR(as.data.frame(as.matrix(mice_merged[["RNA"]]$data)), 
                          ref_mice, ref_mice$label.main)
mice_merged$SingleR.labels <- SingleR_result$labels
DimPlot(mice_merged, reduction = "umap", label = TRUE, pt.size = 0.5, 
        group.by = "SingleR.labels")

################################
# plot cluster proportion
pt <- table(mice_merged$seurat_clusters, mice_merged$orig.ident)
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
Idents(mice_merged) <- "seurat_clusters"
T_cell.markers <- FindConservedMarkers(mice_merged, ident.1 = "12", 
                                       grouping.var = "orig.ident")
mice_merged$seurat_clusters <- factor(mice_merged$seurat_clusters, 
                                      levels = levels(mice_merged))

markers.to.plot <- c("S100a11")
DotPlot(mice_merged, features = markers.to.plot, cols = c("blue", "red"), 
        dot.scale = 8, split.by = "orig.ident") + RotatedAxis()

################################
# Identify DEG across conditions
DEGs_result <- DEGs_across_condition(mice_merged, "mice_control")

# DEG visulization
# Extract top 5 DEGs of each cluster, sorted by increasing P-val
genes.to.label <- extract_top_DEGs(DEGs_result, 
                                   cluster_ident = 0, 
                                   num_gene = 5, 
                                   col_name = "p_val", 
                                   decreasing = F, 
                                   "mice_treatment")
plots <- VlnPlot(mice_merged, features = genes.to.label, 
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
mice_control <- subset(mice_merged, subset = orig.ident == "mice_control")
# Test DEGs identified from Wang_et_al Yuzwa data

# Microglia
FeaturePlot(mice_control, features = Wang_et_al_Marker$sup_table3b$Microglia, 
            reduction = "umap") + 
  plot_annotation(title = "Microglia")

# Astrocyte
FeaturePlot(mice_control, features = Wang_et_al_Marker$sup_table3b$Astro_df, 
            reduction = "umap") + 
  plot_annotation(title = "Astrocyte_common")

# OPC
FeaturePlot(mice_control, features = Wang_et_al_Marker$sup_table3b$OPC, 
            reduction = "umap") + 
  plot_annotation(title = "OPC")

# Oligo
FeaturePlot(mice_control, features = Wang_et_al_Marker$sup_table3b$Oligodendrocyte, 
            reduction = "umap") + 
  plot_annotation(title = "Oligodendrocyte")

# Pre-Oligo
FeaturePlot(mice_control, features = Wang_et_al_Marker$sup_table3b$`Oligodendrocyte-Immature`, 
            reduction = "umap") + 
  plot_annotation(title = "Pre-oligo")

# Vascular
FeaturePlot(mice_control, features = Wang_et_al_Marker$sup_table3b$Vascular, 
            reduction = "umap") + 
  plot_annotation(title = "Vascular")

# RG
FeaturePlot(mice_control, features = Wang_et_al_Marker$sup_table3b$RG_df, 
            reduction = "umap") + 
  plot_annotation(title = "RG-Common")

# Unknown
FeaturePlot(mice_control, features = Wang_et_al_Marker$sup_table3b$Unknown, 
            reduction = "umap") + 
  plot_annotation(title = "Unknown")

# Tri_IPC
FeaturePlot(mice_control, features = Wang_et_al_Marker$sup_table3b$`IPC-Glia`, 
            reduction = "umap") + 
  plot_annotation(title = "Tri_IPC")

# CR
FeaturePlot(mice_control, features = Wang_et_al_Marker$sup_table3b$`Cajal-Retzius cell`, 
            reduction = "umap") + 
  plot_annotation(title = "Cajal-Retzius cell")

# EN
FeaturePlot(mice_control, features = Wang_et_al_Marker$sup_table3b$EN_df, 
            reduction = "umap") + 
  plot_annotation(title = "EN")

# IN
FeaturePlot(mice_control, features = Wang_et_al_Marker$sup_table3b$IN_df, 
            reduction = "umap") + 
  plot_annotation(title = "IN")
################################

# Export mice_merged seurat obj into CSV file
data_to_write_out <- as.data.frame(as.matrix(mice_merged[["RNA"]]$data))
fwrite(x = data_to_write_out, row.names = TRUE, file = "mice_merged.csv")