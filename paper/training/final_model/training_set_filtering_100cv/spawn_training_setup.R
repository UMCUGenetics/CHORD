library(mltoolkit)
library(mutSigExtractor)

options(stringsAsFactors=F)

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

#========= Load data =========#
annotation <- read.delim(paste0(base_dir,'/HMF_DR010_DR047/scripts/det_brca_status/hmf_brca_response.txt'))
contexts <- readRDS(paste0(base_dir,'/HMF_DR010_DR047/matrices/merged_contexts.rds'))

#========== Prep full dataset =========#
## Contexts
full_set <- (function(){
  df <- transformContexts(
    contexts,
    simplify.types=c('snv','indel'),
    rel.types=c('snv','indel','sv')
  )

  ## Add response
  df$response <- annotation$response[match(rownames(df), annotation$sample)] %>% as.factor()
  return(df)
})()

#========= Create training sets =========#
training_samples <- annotation[annotation$in_training_set,'sample']
training_set <- full_set[rownames(full_set) %in% training_samples,]

#========= Create training subsets and holdout 'full' set =========#
spawnTrainingSetups <- function(
  training_set, test_set=NULL,
  seeds=1,
  parent.dir=paste0(base_dir,'/HMF_DR010_DR047/scripts/chord_training/assess_training_set/output/'),
  train.chord.script.path=paste0(base_dir,'/HMF_DR010_DR047/scripts/chord_training/assess_training_set/scripts/train_chord_template.R')
){
  
  ## Create training set
  write.table(training_set, paste0(parent.dir,'/training_set.txt'), sep='\t', quote=F)
  if(!is.null(test_set)){ write.table(test_set, paste0(parent.dir,'/test_set.txt'), sep='\t', quote=F) }
  write.table(annotation, paste0(parent.dir,'/annotation.txt'), sep='\t', quote=F)
  
  ## Create output dir and custom scripts
  pb <- txtProgressBar(min=0, max=length(seeds), initial=0, style=3, width=100)
  counter <- 0
  
  for(i in seeds){
    counter <- counter+1
    setTxtProgressBar(pb,counter)
    
    #i=1
    ## Make dir
    seed_dir <- paste0(parent.dir,'/replicates/seed_',formatC(i, width=nchar(max(seeds)), format="d", flag="0"),'/')
    if(!dir.exists(seed_dir)){ dir.create(seed_dir,recursive=T) }
    
    script_out_path <- paste0(seed_dir,'/train_chord_custom.R')
    file.copy(train.chord.script.path, script_out_path, overwrite=T)
    
    system(sprintf(
      "sed -i '' 's|%s|%s|' %s",
      'PARENT_DIR', 
      sub('/Users/lnguyen','',parent.dir), 
      script_out_path
    ))
    
    system(sprintf(
      "sed -i '' 's|%s|%s|' %s",
      'SEED_DIR', 
      sub('/Users/lnguyen','',seed_dir), 
      script_out_path
    ))
    
    system(sprintf(
      "sed -i '' 's|%s|%s|' %s",
      'SEED', i, script_out_path
    ))
  }
}

#========= Exec =========#
spawnTrainingSetups(training_set, seeds=1:100, parent.dir=paste0(base_dir,'/HMF_DR010_DR047/scripts/chord_training/main/output/60_snvContext_svContext_allSigsRel_noSuppEvidence_customBlacklist3_noBoruta/cv100/'))
