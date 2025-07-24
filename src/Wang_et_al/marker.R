library(readxl)
# This file contain markers utilized in Wang_et_al
# (https://doi.org/10.1038/s41586-024-08351-7)
Wang_et_al_Marker <- list()

##################################################
# Top 4 marker gene documented in Sup table 3b
df <- read_excel("data/raw/Wang_et_al/Supplementary Table 3.xlsx", 
                 sheet = 2, 
                 col_names = F)
colnames(df) <- as.character(df[2, ])
df <- df[-c(1,2), ]

df <- split(df, df$`Cell type`)

# Merge Astrocyte marker from three astrocyte cluster
Astro_df <- merge.data.frame(df$`Astrocyte-Fibrous`, df$`Astrocyte-Immature`, all = T)
Astro_df <- merge.data.frame(Astro_df, df$`Astrocyte-Protoplasmic`, all = T)

# Merge RG marker from three RG cluster
RG_df <- merge.data.frame(df$`RG-oRG`, df$`RG-tRG`, all = T)
RG_df <- merge.data.frame(RG_df, df$`RG-vRG`, all = T)

# Merge EN marker from EN clusters
EN_df <- merge.data.frame(df$`EN-IT-Immature`, df$`EN-L2_3-IT`, all = T)
EN_df <- merge.data.frame(EN_df, df$`EN-L4-IT`, all = T)
EN_df <- merge.data.frame(EN_df, df$`EN-L5-ET`, all = T)
EN_df <- merge.data.frame(EN_df, df$`EN-L5-IT`, all = T)
EN_df <- merge.data.frame(EN_df, df$`EN-L5_6-NP`, all = T)
EN_df <- merge.data.frame(EN_df, df$`EN-L6-CT`, all = T)
EN_df <- merge.data.frame(EN_df, df$`EN-L6-IT`, all = T)
EN_df <- merge.data.frame(EN_df, df$`EN-L6b`, all = T)
EN_df <- merge.data.frame(EN_df, df$`EN-Newborn`, all = T)
EN_df <- merge.data.frame(EN_df, df$`EN-Non-IT-Immature`, all = T)
EN_df <- merge.data.frame(EN_df, df$`IPC-EN`, all = T)

# Merge IN marker from IN clusters
IN_df <- merge.data.frame(df$`IN-CGE-Immature`, df$`IN-CGE-SNCG`, all = T)
IN_df <- merge.data.frame(IN_df, df$`IN-CGE-VIP`, all = T)
IN_df <- merge.data.frame(IN_df, df$`IN-dLGE-Immature`, all = T)
IN_df <- merge.data.frame(IN_df, df$`IN-MGE-Immature`, all = T)
IN_df <- merge.data.frame(IN_df, df$`IN-MGE-PV`, all = T)
IN_df <- merge.data.frame(IN_df, df$`IN-MGE-SST`, all = T)
IN_df <- merge.data.frame(IN_df, df$`IN-Mix-LAMP5`, all = T)
###############################
# Merge duplicated gene
for (df2 in c("Astro_df", "RG_df", "EN_df", "IN_df")) {
  # Step 1: Find duplicated Genes (i.e., Genes with multiple entries)
  dup_genes <- get(df2)$Gene[duplicated(get(df2)$Gene) | duplicated(get(df2)$Gene, fromLast = TRUE)]
  
  # Step 2: Merge duplicates
  merged_dup_df <- get(df2) %>%
    filter(Gene %in% dup_genes) %>%
    mutate(across(-Gene, ~ as.numeric(as.character(.)), .names = "{.col}")) %>%
    group_by(Gene) %>%
    summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop") %>%
    mutate(duplicated = TRUE) %>%
    select(duplicated, Gene, everything())
  merged_dup_df <- merged_dup_df[,-3]
  
  # Step 3: Retain non-duplicated rows
  nondup_df <- get(df2) %>%
    filter(!(Gene %in% dup_genes)) %>%
    mutate(across(-Gene, ~ as.numeric(as.character(.)), .names = "{.col}")) 
  
  # Step 4: Combine both
  final_df <- bind_rows(merged_dup_df, nondup_df)
  df[[df2]] <- final_df
}
rm(EN_df)
rm(IN_df)
rm(RG_df)
rm(Astro_df)
rm(final_df)
rm(merged_dup_df)
rm(nondup_df)
rm(df2)
rm(dup_genes)
####################################################
# trimmed_list <- lapply(df, function(df) head(df[, 2, drop = FALSE], 4))
trimmed_list <- lapply(df, function(d) {
  d$Average_log2_fold_change <- as.numeric(as.character(d$Average_log2_fold_change))
  d <- d[order(-d$Average_log2_fold_change), ]
  d[1:4, 2, drop = FALSE]
})

#entry_list <- lapply(df, function(df) {
#  head(as.character(df[[2]]), 4)
#})

entry_list <- lapply(df, function(d) {
  d$Average_log2_fold_change <- as.numeric(as.character(d$Average_log2_fold_change))
  d <- d[order(-d$Average_log2_fold_change), ]
  head(as.character(d[[2]]), 4)
})

capitalize_first_list <- function(entry_list) {
  lapply(entry_list, function(vec) {
    sapply(vec, function(str) {
      str <- tolower(str)
      if (nchar(str) == 0) return(str)
      paste0(toupper(substr(str, 1, 1)), substr(str, 2, nchar(str)))
    }, USE.NAMES = FALSE)
  })
}

capitalized <- capitalize_first_list(entry_list)

Wang_et_al_Marker[["sup_table3b"]] <- capitalized

rm(df)
rm(trimmed_list)
rm(entry_list)
rm(capitalized)
##################################################
# marker gene from snMultiome analysis shown on paper
snMultiome_marker <- list()

snMultiome_marker[["Microglia"]] <- c("Irf8")

snMultiome_marker[["oligo"]] <- c("Grm3", "Olig1", "Sox10", "Mbp")

snMultiome_marker[["OPC"]] <- c("Map3k1", "Olig1", "Sox10", "Bcas1", "Mbp")

snMultiome_marker[["Astro"]] <- c("Gfap", "Aqp4", "Moxd1", "Grm3")

snMultiome_marker[["RG"]] <- c("Hes1", "Tfap2c", "Cryab", "Fbxo32", "Tnc", "Hopx")

snMultiome_marker[["EN"]] <- c("Slc17a6", "Eomes", "Nrp1", "Glis3", "Cux2", "Rorb", 
               "Il1rapl2", "Znf804b", "Fezf2", "Nwd2", "Tshz2", "Syt6", 
               "Fam160a1")
snMultiome_marker[["IN"]] <- c("Gad1", "Meis2", "Adarb2", "Vip", "Cck", "Lamp5", "Nxph1", 
               "Sst", "Pvalb")
snMultiome_marker[["Vascular"]] <- c("Igfbp7")

snMultiome_marker[["Tri_IPC"]] <- c("Egfr", "Map3k1", "Olig1")

snMultiome_marker[["CR"]] <- c("Reln")

Wang_et_al_Marker[["snMultiome_marker"]] <- snMultiome_marker

rm(snMultiome_marker)
################################################
# Surface protein utilized for FACS
FACS_marker <- list()

FACS_marker[["Tri_IPC"]] <- c("Egfr", "F3", "Pdgfra")

FACS_marker[["RG"]] <- c("Itga2", "Egfr")

FACS_marker[["OPC"]] <- c("Pdgfra")

Wang_et_al_Marker[["FACS_marker"]] <- FACS_marker

rm(FACS_marker)
################################################
# Marker utilized at DIV1 in vitro cell culture
Vitro_marker <- list()

Vitro_marker[["RG"]] <- c("Tfap2c", "Cryab")

Vitro_marker[["Tri_IPC"]] <- c("Olig2", "Egfr")

Wang_et_al_Marker[["Vitro_marker"]] <- Vitro_marker

rm(Vitro_marker)