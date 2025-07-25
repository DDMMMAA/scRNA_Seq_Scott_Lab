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
    message("#   start: the start of range of resolution (>= 0, default == 0)")
    message("#   end: the end of range of resolution (>= start, default == 1.5)")
    message("#   margin: the margin of range (positive rational number, default == 0.1)")
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
    message("#   num_top_gene: number of top differentially expressed gene used for GO (>0, default == 30)")
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
# wrapper for DEG of clusters across conditions
# Perform DEG of "control" of all clusters w/ respect to all treatment group
# Return a list of list of data.frame that store DEGs
# Argument:
#   obj: A seurat object
#   control_name: the name of "Control" group used in orig.ident
DEGs_across_condition <- function(obj, control_name) {
  if(any(class(obj) != "Seurat") ||
    (class(control_name) != "character")){
    message("# Argument:")
    message("#   obj: a Seurat object")
    message("#   control_name: the name of 'Control' group used in orig.ident")
  } else {
    obj$celltype.condition <- paste(obj$seurat_clusters, 
                                            obj$orig.ident, sep = "_")
    Idents(obj) <- "celltype.condition"
    result <- list()
    treatment_groups <- unique(obj$orig.ident)
    treatment_groups <- treatment_groups[ !treatment_groups == control_name]
    for (i in levels(obj$seurat_clusters)) {
      # The line above might be redundant, comment out to test
      # cluster_name <- paste("Cluster", i, sep = "")
      result_name <- paste("DEG_result_cluster", i, sep = "")
      enrich_res <- list()
      message(treatment_groups)
      for (treatment_group in treatment_groups) {
        enrich_res[[treatment_group]] <- FindMarkers(obj, 
                                  ident.1 = paste(i, control_name, sep = "_"), 
                                  ident.2 = paste(i, treatment_group, sep = "_"))
        # This line of code calculate the expected value of FC value is true
        # Under the Ha assumption for each gene. It doesn't really have a biological
        # interpretation (at least I can't think about one), but an two variable 
        # function of p-value and avg_log2FC that generate a value to rank the 'significance'
        # of DEGs. The propose of this value is 
        # 1: disregard the ambiguous p-value == 0.05 absolute threshold 
        # (https://doi.org/10.1080/00031305.2016.1154108)
        # 2: Reduce the subjective human bias of determine the
        # 'significance' of DEGs
        # 3: Save time so that I don't need to worry p-value and avg_log2FC
        # respectively
        enrich_res[[treatment_group]]$p_FC <- ((1 - enrich_res[[treatment_group]]$p_val) * 
                                                 (enrich_res[[treatment_group]]$avg_log2FC))
      }
      result[[result_name]] <- enrich_res
      message("Complete cluster: ", i)
    }
    return(result)
  }
}

###################################
# Wrapper to extract top DEGs
# rank top DEGs based on either P-val or avg_log2FC
# Return a vector of gene symbol
# Argument:
#   obj: a list of list of dataframe, each dataframe contain gene symbol as row names, 
#        p-val & avg_log2FC as column
#   cluster_ident: the numeric identity of a cluster (cluster_ident < length(obj))
#   num_gene: the number of top DEGs to extract (default == 5)
#   col_name: either "p_val" or "avg_log2FC" ("p_val" || "avg_log2FC")
#   decreasing: order of col_name, decreasing or not (T/F, default == F)
#   treatment_group: The group of DEGs w/ respect to
extract_top_DEGs <- function(obj, 
                             cluster_ident, 
                             num_gene = 5, 
                             col_name, 
                             decreasing = F, 
                             treatment_group) {
  if ((class(obj) != "list") || 
      ((class(cluster_ident) != "numeric") && (class(cluster_ident) != "integer")) ||
      (cluster_ident >= length(obj)) || 
      ((class(num_gene) != "numeric") && (class(num_gene) != "integer")) ||
      (!col_name %in% c("p_val", "avg_log2FC", "p_FC")) ||
      (class(decreasing) != "logical") ||
      (class(treatment_group) != "character")) {
    message("# Argument:")
    message("#   obj: a list of dataframe, each dataframe contain gene symbol as row names, ")
    message("#        p-val & avg_log2FC as column")
    message("#   cluster_ident: the numeric identity of a cluster (cluster_ident < length(obj))")
    message("#   num_gene: the number of top DEGs to extract (default == 5)")
    message("#   col_name: either 'p_val' or 'avg_log2FC' ('p_val' || 'avg_log2FC')")
    message("#   decreasing: order of col_name, decreasing or not (T/F, default == F)")
  } else {
    cluster_str <- paste("DEG_result_cluster", cluster_ident, sep = "")
    return(head(row.names(obj[[cluster_str]][[treatment_group]]
                          [order(obj[[cluster_str]][[treatment_group]][[col_name]], 
                                 decreasing = decreasing), ]), num_gene))
  }
}

################################
# Wrapper for subsetting DEGs
# Return a list of list of data.frame that store DEGs
# Argument:
#   obj: a list of list of data.frame that store DEGs
#   logic_exp: a logical expression that determine which DEGs to subset ("character")
subset_DEGs <- function(obj, logic_exp) {
  if ((class(obj) != "list") || 
      (class(logic_exp) != "character")) {
    message("# Argument:")
    message("#   obj: a list of list of dataframe storing DEGs")
    message("#   logic_exp: a logical expression that determine which DEGs to subset ('character')")
  } else {
    result_cluster <- list()
    for (i in names(obj)) {
      result_cluster[[i]] <- list()
      for (j in names(obj[[i]])) {
        df <- obj[[i]][[j]]
        result_cluster[[i]][[j]] <- df[eval(parse(text = logic_exp), envir = df), ]
      }      
    }
    return(result_cluster)
  }
}

################################
# Wrapper for run enrichGO among all DEGs, export GO barplot optionally
# Return a list of list of data.frame that store DEGs' GO result
# Argument:
#   obj: a list of list of data.frame that store DEGs
#   export_figure: export GO barplot figure or not (T/F)
#   dir: the dir of exported figure, if "export_figure" == T ("character")
#   width: width of exported figure, if "export_figure" == T ("numeric")
#   height: height of exported figure, if "export_figure" == T ("numeric")
DEGs_enrichGO <- function(obj, export_figure = F, dir, width = 1280, height = 720) {
  if ((class(obj) != "list") || 
      (class(export_figure) != "logical") ||
      (class(dir) != "character") ||
      (class(width) != "numeric") ||
      (class(height) != "numeric")) {
    message("# Argument:")
    message("#   obj: a list of list of data.frame that store DEGs")
    message("#   export_figure: export GO barplot figure or not (T/F)")
    message("#   dir: the dir of exported figure, if 'export_figure' == T ('character')")
    message("#   width: width of exported figure, if 'export_figure' == T ('numeric')")
    message("#   height: height of exported figure, if 'export_figure' == T ('numeric')")
  } else {
    GO_result_cluster <- list()
    # Test is figure need to be export and provided dir doesn't exists
    if ((export_figure == T) && (!dir.exists(dir))) {
      dir.create(dir)
    }
    for (i in names(obj)) {
      GO_result_cluster[[i]] <- list()
      # Test is figure need to be export and cluster_subdir of provided dir doesn't exists
      cluster_subdir <- paste(dir, i, sep = "/")
      if ((export_figure == T) && (!dir.exists(cluster_subdir))) {
        dir.create(cluster_subdir)
      }
      for (j in names(obj[[i]])) {
        df <- obj[[i]][[j]]
        genes.to.GO <- row.names(df)
        genes.to.GO <- mapIds(org.Mm.eg.db, keys = genes.to.GO, keytype = "SYMBOL", column="ENSEMBL")
        GO_result_cluster[[i]][[j]] <- enrichGO(gene = genes.to.GO, 
                                                keyType = "ENSEMBL", 
                                                OrgDb = "org.Mm.eg.db", 
                                                ont = "BP")
        plot <- barplot(GO_result_cluster[[i]][[j]], showCategory = 10)
        treatment_GO_name <- paste(cluster_subdir, j, sep = "/")
        treatment_GO_name <- paste(treatment_GO_name, ".png", sep = "")
        # Test if figure need to be export
        if((export_figure == T) && 
           (nrow(subset(GO_result_cluster[[i]][[j]]@result, subset = (p.adjust <= 0.05)))) != 0) {
          # export figure
          png(filename = treatment_GO_name, width = width, height = height)
          print(plot)
          dev.off()
        }
        message("Complete treatment:", j)
      }
      message("Complete Cluster:", i)
    }
    return(GO_result_cluster)
  }
}

#######################################
# Wrapper to merge two DEG data.frames by gene (row name), averaging all numeric columns
# Return a merged data.frame with gene names as rownames and averaged columns
# Argument:
#   df1: first DEG data.frame (row names = gene symbols; all columns numeric)
#   df2: second DEG data.frame (same structure as df1)
#   keep_non_common: whether to keep genes that are not shared (TRUE/FALSE)
merge_DEG_dfs <- function(df1, df2, keep_non_common = FALSE) {
  # Error handling
  if (!is.data.frame(df1) || !is.data.frame(df2)) {
    message("# Argument:")
    message("#   df1: first DEG data.frame (row names = gene symbols; all columns numeric)")
    message("#   df2: second DEG data.frame (same structure as df1)")
    message("#   keep_non_common: whether to keep genes that are not shared (TRUE/FALSE)")
    stop("Both df1 and df2 must be data.frames.")
  }
  
  if (!identical(colnames(df1), colnames(df2))) {
    stop("The two data.frames do not have the same column names.")
  }
  
  # Ensure all values are numeric
  df1[] <- lapply(df1, function(x) as.numeric(as.character(x)))
  df2[] <- lapply(df2, function(x) as.numeric(as.character(x)))
  
  # Identify genes
  genes_df1 <- rownames(df1)
  genes_df2 <- rownames(df2)
  common_genes <- intersect(genes_df1, genes_df2)
  
  # Merge common genes
  df1_common <- df1[common_genes, , drop = FALSE]
  df2_common <- df2[common_genes, , drop = FALSE]
  merged_common <- (df1_common + df2_common) / 2
  rownames(merged_common) <- common_genes
  
  # Handle non-common genes
  if (keep_non_common) {
    unique_df1 <- df1[setdiff(genes_df1, common_genes), , drop = FALSE]
    unique_df2 <- df2[setdiff(genes_df2, common_genes), , drop = FALSE]
    merged_result <- rbind(merged_common, unique_df1, unique_df2)
  } else {
    merged_result <- merged_common
  }
  
  message("Merged ", length(common_genes), " shared genes.",
          if (keep_non_common) paste0(" Included ", 
                                      nrow(merged_result) - length(common_genes), 
                                      " non-common genes."))
  return(merged_result)
}
