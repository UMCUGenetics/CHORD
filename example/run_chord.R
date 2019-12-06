library(CHORD)

base_dir <- '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/CHORD/example/'

options(stringsAsFactors=F)

#========= Load vcfs =========#
## Get vcf file paths
vcf_files <- data.frame(
  snv_indel=list.files(paste0(base_dir,'/vcf/'), pattern='*snv_indel.vcf.gz', full.names=T),
  sv=list.files(paste0(base_dir,'/vcf/'), pattern='*sv.vcf.gz', full.names=T)
)

## Assign sample names
vcf_files <- cbind(
  sample=sapply(strsplit(basename(vcf_files$snv_indel),'_'),`[`,1),
  vcf_files
)

#========= Extract mutation contexts =========#
## Make dir to output contexts
contexts_dir <- paste0(base_dir,'/contexts/')
dir.create(contexts_dir)
dir.create(paste0(contexts_dir,'/raw/'))

## Extract contexts for all samples
for(i in 1:nrow(vcf_files)){
  #i=1
  params <- as.list(vcf_files[i,])
  out_path <- paste0(contexts_dir,'/raw/',params$sample,'_contexts.txt')
  
  if(!file.exists(out_path)){
    
    message('Extracting contexts for: ',params$sample)
    
    extractSigsChord(
      vcf.snv=params$snv_indel,
      vcf.sv=params$sv, sv.caller='manta',
      sample.name=params$sample,
      output.path=out_path, verbose=F
    )
    
  } else {
    message('Contexts exist for: ',params$sample, ' (Skipping extracting contexts)')
  }
  
}

## Merge contexts into one dataframe
merged_contexts_path <- paste0(contexts_dir,'/merged_contexts.txt')
if(!file.exists(merged_contexts_path)){
  context_files <- list.files(paste0(contexts_dir,'/raw/'), full.names=T)
  merged_contexts <- do.call(rbind, lapply(context_files, function(i){
    read.delim(i, check.names=F)
  }))
  write.table(merged_contexts, merged_contexts_path, sep='\t', quote=F)
} else {
  merged_contexts <- read.delim(merged_contexts_path, check.names=F)
}

#========= Predict HRD with CHORD =========#
pred <- chordPredict(merged_contexts, verbose=F)









