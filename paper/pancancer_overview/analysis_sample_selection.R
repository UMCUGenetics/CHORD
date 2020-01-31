#========= Path prefixes =========#
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
metadata <- read.delim(paste0(base_dir,'/HMF_DR010_DR047/metadata/DR-047-update1_metadata.txt'))

#--------- HMF predictions ---------#
## Features
contexts_hmf <- readRDS(paste0(base_dir,'/HMF_DR010_DR047/matrices/merged_contexts.rds'))

#========= Selection based on max purity =========#
multi_biop_patients <- metadata$patient_id[ duplicated(metadata$patient_id) ]

metadata_ss <- metadata[!(metadata$patient_id %in% multi_biop_patients),]

metadata_ss_multi_biop <- do.call(rbind,lapply(multi_biop_patients, function(i){
   df <- metadata[metadata$patient_id==i,]
   df[which.max(df$tumor_purity),]
}))

metadata_ss <- rbind(metadata_ss, metadata_ss_multi_biop)

df <- metadata[c('sample_id','hmf_id','patient_id','primary_tumor_location','tumor_purity')]
df$is_max_purity_biopsy <- df$sample_id %in% metadata_ss$sample_id

#========= Selection based on mut load =========#
indel_load <- rowSums(contexts_hmf$indel)

df$indel_load <- indel_load[ match(df$sample_id, names(indel_load)) ]
df$indel_load_ge50 <- df$indel_load >= 50

sv_load <- rowSums(contexts_hmf$sv)
#df$sv_load <- indel_load[ match(df$sample_id, names(sv_load)) ]
df$sv_load <- sv_load[ match(df$sample_id, names(sv_load)) ]
df$sv_load_ge30 <- df$sv_load >= 30


indel.rep_load <- rowSums(contexts_hmf$indel[,grep('rep',colnames(contexts_hmf$indel))])
df$indel.rep_load <- indel.rep_load[ match(df$sample_id, names(indel.rep_load)) ]
df$has_msi <- df$indel.rep_load > 14000


#========= Final selection + export =========#
## Main
df$is_selected <- df$indel_load_ge50 & df$sv_load_ge30 & df$is_max_purity_biopsy & !df$has_msi

df <- df[
   with(df,{ order(patient_id, tumor_purity) })
,]

write.table(
   df,
   paste0(base_dir,'/HMF_DR010_DR047/analysis/pancancer_overview/data/analysis_sample_selection.txt'),
   sep='\t',row.names=F,quote=F
)


## Supplementary table
df_supp_table <- df

df_supp_table$patient_id <- gsub('\\w$','',df_supp_table$hmf_id)
df_supp_table[c('sample_id','primary_tumor_location')] <- NULL

colnames(df_supp_table)[1] <- 'sample_id'

write.table(
   df_supp_table,
   paste0(base_dir,'/HMF_DR010_DR047/analysis/pancancer_overview/supp_tables/pancancer_analysis_sample_selection.txt'),
   sep='\t',row.names=F,quote=F
)
