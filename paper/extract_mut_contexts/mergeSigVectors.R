library(mutSigExtractor)

base_dir <- '/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/'
if(dir.exists('/Users/lnguyen/')){
   base_dir <- paste0('/Users/lnguyen/', base_dir)
}

## Load contexts
context_dirs <- list(
   snv = paste0(base_dir,'/HMF_DR010_DR047/matrices/snv/'),
   indel = paste0(base_dir,'/HMF_DR010_DR047/matrices/indel/'),
   sv = paste0(base_dir,'/HMF_DR010_DR047/matrices/sv/')
)
contexts <- lapply(context_dirs, readSigsAsDf)

## Get common samples (for when any snv/indel/sv files are missing)
common_samples <- Reduce(intersect, lapply(contexts,rownames))

contexts_ss <- lapply(contexts, function(i){
  i[rownames(i) %in% common_samples,]
})

saveRDS(contexts_ss, paste0(base_dir,'/HMF_DR010_DR047/matrices/merged_contexts.rds'))
write.table(do.call(cbind, unname(contexts_ss)), paste0(base_dir,'/HMF_DR010_DR047/matrices/hmf_contexts.txt'), sep='\t', quote=F)
