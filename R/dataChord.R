#' CHORD (Classifier for Homologous Recombination Deficiency)
#'
#' A random forest predictor of homologous recombination deficiency
#'
#' Version 60: 60_snvContext_svContext_allSigsRel_noSuppEvidence_customBlacklist3_noBoruta
#'
#' @docType data
#'
#' @usage data(CHORD)
'CHORD'

# #--------- Code used to save the random forest model into the CHORD R package ---------#
# model_dir <- '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_DR010_DR047/scripts/chord_training/main/output/60_snvContext_svContext_allSigsRel_noSuppEvidence_customBlacklist3_noBoruta/seed015/'
# CHORD <- readRDS(paste0(model_dir,'/final/rf_out.rds'))$RF
# save(CHORD,file='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/CHORD/data/CHORD.rda')
