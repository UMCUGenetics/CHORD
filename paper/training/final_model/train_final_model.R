#!/usr/bin/env Rscript

options(stringsAsFactors=F)

#========= Packages =========#
if(.Platform$GUI == 'RStudio'){
  library(mltoolkit)
  library(mutSigExtractor)
} else {
  library(devtools)
  load_all('/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/mltoolkit/')
  load_all('/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/mutSigExtractor/')
  load_all('/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/CHORD/')
}

library(Boruta)
library(doMC)

#========= Paths =========#
base_dir <- '/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/'
if(dir.exists('/Users/lnguyen/')){
  base_dir <- paste0('/Users/lnguyen/', base_dir)
}

contexts <- list(
  hmf=readRDS(paste0(base_dir,'/HMF_DR010_DR047/matrices/merged_contexts.rds')),
  brca_eu=readRDS(paste0(base_dir,'/BRCA_EU/matrices/main/merged_contexts.rds'))
)

annotation <- read.delim(paste0(base_dir,'/HMF_DR010_DR047/scripts/det_brca_status/hmf_brca_response.txt'))

#========= Make signature matrix =========#
#--------- SNV,SV contexts; all rel ---------#
out_dir <- paste0(base_dir,'/HMF_DR010_DR047/scripts/chord_training/main/output/60_snvContext_svContext_allSigsRel_noSuppEvidence_customBlacklist3_noBoruta/')
prepInput <- function(contexts){
  transformContexts(
    contexts,
    simplify.types=c('snv','indel'),
    rel.types=c('snv','indel','sv')
  )
}

#--------- HMF ---------#
hmf_full_set <- (function(){
  df <- as.data.frame(prepInput(contexts$hmf))
  df$response <- annotation$response[match(rownames(df), annotation$sample)] %>% as.factor()
  return(df)
})()

## Blacklisting
custom_blacklist <- unlist(readRDS(paste0(base_dir,'/HMF_DR010_DR047/scripts/chord_training/main/output/60_snvContext_svContext_allSigsRel_noSuppEvidence_customBlacklist3_noBoruta/cv100/custom_blacklist.rds')))
annotation[annotation$sample %in% custom_blacklist,'in_training_set'] <- FALSE

## Make 
training_samples <- annotation[annotation$in_training_set,'sample']
hmf_training_set <- hmf_full_set[rownames(hmf_full_set) %in% training_samples,]

#--------- Export matrices ---------#
write.tsv <- function(data, ...){ write.table(data, sep='\t', quote=F, ...) }
write.tsv(hmf_full_set, paste0(out_dir,'/hmf_full_set.txt'))
write.tsv(hmf_training_set, paste0(out_dir,'/hmf_training_set.txt'))
write.tsv(annotation, paste0(out_dir,'/hmf_annotation.txt'), row.names=F)

#========= Execute =========#
#--------- Save R script ---------# #!! Only works when current file is sourced
file.copy(
  sys.frame(1)$ofile,
  to = file.path(out_dir, paste0('train_chord_', Sys.Date(), ".R"))
)

par <- list(
  seed = 15,
  do.feature.prefilter = T,
  feature.prefilter.p.value.n = 1,
  feature.prefilter.p.value.p = 0.01,
  do.balance = T, 
  do.boruta = F
)

#--------- Final model ---------#
randomForestTrainTest.chord(
  train = hmf_training_set,
  export.dir = paste0(out_dir,'/final/'),
  train.test.seed = par$seed,
  do.feature.prefilter = par$do.feature.prefilter,
  feature.prefilter.p.value.n = par$feature.prefilter.p.value.n,
  feature.prefilter.p.value.p = par$feature.prefilter.p.value.p,
  do.balance = par$do.balance,
  do.boruta = par$do.boruta
)


