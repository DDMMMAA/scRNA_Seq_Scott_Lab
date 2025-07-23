# This file contain markers utilized in Wang_et_al
# (https://doi.org/10.1038/s41586-024-08351-7)
Wang_et_al_Marker <- list()

##################################################
# marker identified via snMultiome analysis
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

Vitro_marker[["Tri_IPC"]] <- c("Olig2", "EGFR")

Wang_et_al_Marker[["Vitro_marker"]] <- Vitro_marker

rm(Vitro_marker)