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

options(stringsAsFactors=F)

#========= Paths =========#
base_dir <- '/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/'
if(dir.exists('/Users/lnguyen/')){
  base_dir <- paste0('/Users/lnguyen/', base_dir)
}

parent_dir <- 'PARENT_DIR'
seed_dir <- 'SEED_DIR'

#========= Load data =========#
training_set <- read.delim(paste0(parent_dir,'/training_set.txt'), stringsAsFactors=T)

#========= Train + predict =========#
par <- list(
  seed = SEED,
  do.feature.prefilter = T,
  feature.prefilter.p.value.n = 1,
  feature.prefilter.p.value.p = 0.01,
  do.balance = T,
  do.boruta = F
)

randomForestCvNested.chord(
  df = training_set,
  parent.export.dir = paste0(seed_dir,'/cv/'),
  rfcv.nested.seed = par$seed,
  do.feature.prefilter = par$do.feature.prefilter,
  feature.prefilter.p.value.n = par$feature.prefilter.p.value.n,
  feature.prefilter.p.value.p = par$feature.prefilter.p.value.p,
  do.balance = par$do.balance,
  do.boruta = par$do.boruta
)
