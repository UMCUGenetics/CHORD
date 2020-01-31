#!/usr/bin/env Rscript

options(stringsAsFactors=F)

#========= Packages =========#
if(.Platform$GUI == 'RStudio'){
  library(mltoolkit)
  library(mutSigExtractor)
  library(CHORD)
} else {
  library(devtools)
  load_all('/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/mltoolkit/')
  load_all('/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/mutSigExtractor/')
  load_all('/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/CHORD/')
}

library(Boruta)
library(doMC)

#========= Paths =========#
base_dir <- list(
  hpc='/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/',
  mnt='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/',
  local='/Users/lnguyen/Documents/Luan_projects/CHORD/'
)

for(i in base_dir){
  if(dir.exists(i)){ 
    base_dir <- i 
    break
  }
}

#========= Make signature matrix =========#
#--------- SNV,SV contexts; all rel ---------#
args <- commandArgs(trailingOnly=T)
out_dir <- args[1]
hmf_training_set <- read.delim(paste0(out_dir,'/training_set.txt'))
hmf_test_set <- read.delim(paste0(out_dir,'/test_set.txt'))

## Blacklisting
custom_blacklist <- as.character(unlist(readRDS(paste0(out_dir,'/custom_blacklist.rds'))))
hmf_training_set <- hmf_training_set[!(rownames(hmf_training_set) %in% custom_blacklist),]

hmf_training_set$response <- as.factor(hmf_training_set$response)

#--------- Export matrices ---------#
write.tsv <- function(data, ...){ write.table(data, sep='\t', quote=F, ...) }

#========= Execute =========#
par <- list(
  seed = 1,
  do.feature.prefilter = T,
  feature.prefilter.p.value.n = 1,
  feature.prefilter.p.value.p = 0.01,
  do.balance = T, 
  do.boruta = F
)

#--------- Final model ---------#
randomForestTrainTest.chord(
  train = hmf_training_set,
  test = hmf_test_set,
  export.dir = paste0(out_dir,'/final/'),
  train.test.seed = par$seed,
  do.feature.prefilter = par$do.feature.prefilter,
  feature.prefilter.p.value.n = par$feature.prefilter.p.value.n,
  feature.prefilter.p.value.p = par$feature.prefilter.p.value.p,
  do.balance = par$do.balance,
  do.boruta = par$do.boruta
)


