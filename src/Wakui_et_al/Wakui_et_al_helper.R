# This file contain script that align barcode in Wakui_et_al's expression matrix
# to treatment group
#############################
# Align orig.ident of given Seurat obj to their corresponding treatment group
# based on provided barcode-treatment group match stored in given dir
# return: a Seurat obj
# Argument:
#   obj: a Seurat object
#   dir: directory that store barcode-treatment group match .tsv file
#   exp_repeat: the repeat number of experiment (1 || 2 || 3 || 4)
Align_barcode_group <- function(obj, dir, exp_repeat) {
  # Precondition checking
  if ((class(obj) != "Seurat") || 
      (class(dir) != "character") ||
      (!exp_repeat %in% c(1, 2, 3, 4))) {
    message("# Argument:")
    message("#   obj: a Seurat object")
    message("#   dir: directory that store barcode-treatment group match .tsv file")
    message("#   exp_repeat: the repeat number of experiment (1 || 2 || 3 || 4)")
  } else {
    CI_barcode_directory <- paste(dir, "/Cortex_CI_", exp_repeat, ".tsv", sep = "")
    HIE_barcode_directory <- paste(dir, "/Cortex_HIE_", exp_repeat, ".tsv", sep = "")
    Hypoxia_barcode_directory <- paste(dir, "/Cortex_Hypoxia_", exp_repeat, ".tsv", sep = "")
    Sham_barcode_directory <- paste(dir, "/Cortex_Sham_", exp_repeat, ".tsv", sep = "")
    
    CI_barcode <- read.table(CI_barcode_directory)
    HIE_barcode <- read.table(HIE_barcode_directory)
    Hypoxia_barcode <- read.table(Hypoxia_barcode_directory)
    Sham_barcode <- read.table(Sham_barcode_directory)
    
    CI_barcode$group <- "CI"
    HIE_barcode$group <- "HIE"
    Hypoxia_barcode$group <- "Hypoxia"
    Sham_barcode$group <- "Sham"
    
    barcode <- merge(CI_barcode, HIE_barcode, all = T)
    barcode <- merge(barcode, Hypoxia_barcode, all = T)
    barcode <- merge(barcode, Sham_barcode, all = T)
    
    message(length(rownames(Wakui_cortex1@meta.data)) == length(barcode$V1))
    message(any(rownames(obj@meta.data) != barcode$V1))
    
    obj[[]]$orig.ident <- barcode$group
    return(obj)
  }
}
