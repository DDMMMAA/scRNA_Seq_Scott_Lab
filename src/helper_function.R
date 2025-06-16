# This file contain all self-defined helper function

######################
# Find cluster of given Seurat obj under range of resolution
# if range is not a integer division of margin, 
# stop at the greatest resolution that less than or equal to end
# return : a Seurat obj
# Argument:
#   obj: a Seurat object
#   start: the start of range of resolution (>= 0)
#   end: the end of range of resolution (>= start)
#   margin: the margin of range (positive rational number)
FindClusters.range <- function(obj, start = 0, end = 1.5, margin = 0.1) {
  # Precondition checking
  if ((class(obj) != "Seurat") || 
      (start < 0) || 
      (end < start) || 
      (margin <= 0)) {
    message("# Argument:")
    message("#   obj: a Seurat object")
    message("#   start: the start of range of resolution (>= 0, default = 0)")
    message("#   end: the end of range of resolution (>= start, default = 1.5)")
    message("#   margin: the margin of range (positive rational number, default = 0.1)")
  } else {
    range <- seq(start, end, margin)
    # loop over range
    for(i in range) {
      obj <- FindClusters(obj, resolution = i)
      message("Added clustering under resolution: ", i)
    }
    return(obj)
  }
}

######################
# Run GO iteratively among clusters of a given data.fram
# Create GO_result vector that store each cluster's GO result at the corresponding index
# Return vector of GOenrichment object
# Argument:
#   obj: A data.frame generated from Seurat FindAllMarkers() function
#   num_top_gene: number of top differentially expressed gene used for GO (>0)
#   ont: One of "BP", "MF", and "CC" subontologies, or "ALL" for all three
enrichGO.clusters <- function(obj, num_top_gene = 30, ont) {
  if ((class(obj) != "data.frame") || 
      (num_top_gene <= 0) || 
      (!ont %in% c("BP", "MF", "CC", "ALL"))) {
    message("# Argument:")
    message("#   obj: a Seurat object")
    message("#   num_top_gene: number of top differentially expressed gene used for GO (>0, default = 30)")
    message("#   ont: One of 'BP', 'MF', and 'CC' subontologies, or 'ALL' for all three")
  } else {
    result <- list()
    for (i in 0:(length(levels(obj$cluster)) - 1)) {
      cluster_name <- paste("Cluster", i, sep = "")
      result_name <- paste("result_cluster", i, sep = "")
      assign(cluster_name, 
             head(subset(obj, cluster == i)[["gene"]], num_top_gene))
      buf <- get(cluster_name)
      enrich_res <- enrichGO(gene = buf, 
                             keyType = "SYMBOL", 
                             OrgDb = "org.Mm.eg.db", 
                             ont = ont)
      result[[result_name]] <- enrich_res
      message("Complete cluster: ", i)
    }
    return(result)
  }
}

##################################
# wrapper for DEG of cluster across conditions
# Perform DEF of all clusters across conditions
# Return a vector of data.frame that store DEGs
# Argument:
#   obj: A seurat object
DEGs_across_condition <- function(obj) {
  if(any(class(obj) != "Seurat")) {
    message("# Argument:")
    message("#   obj: a Seurat object")
  } else {
    obj$celltype.condition <- paste(obj$seurat_clusters, 
                                            obj$orig.ident, sep = "_")
    Idents(obj) <- "celltype.condition"
    result <- list()
    for (i in levels(obj$seurat_clusters)) {
      cluster_name <- paste("Cluster", i, sep = "")
      result_name <- paste("DEG_result_cluster", i, sep = "")
      enrich_res <- FindMarkers(obj, 
                                ident.1 = paste(i, "_mice_control", sep = ""), 
                                ident.2 = paste(i, "_mice_treatment", sep = ""))
      result[[result_name]] <- enrich_res
      message("Complete cluster: ", i)
    }
    return(result)
  }
}
