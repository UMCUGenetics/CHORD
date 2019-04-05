#' CHORD (Classifier for Homologous Recombination Deficiency)
#'
#' A random forest predictor of homologous recombination deficiency
#'
#' Version 38: 38_HMF_pkgMutSigs_germLt3-somEq0-lohEq0_greylistIncl_multiclass_snvContext_svContext_allSigsRel_scaleGridSearch_featPreFilt_boruta_newAnn
#'
#' @docType data
#'
#' @usage data(CHORD)
'CHORD'

#--------- Code used to save the random forest model into the CHORD R package ---------#
# model_dir <- '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/CHORD/training/models/38_HMF_pkgMutSigs_germLt3-somEq0-lohEq0_greylistIncl_multiclass_snvContext_svContext_allSigsRel_scaleGridSearch_featPreFilt_boruta_newAnn/'
# CHORD <- readRDS(paste0(model_dir,'/final/rf_out.rds'))$RF
# save(CHORD,file='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/CHORD/data/CHORD.rda')
