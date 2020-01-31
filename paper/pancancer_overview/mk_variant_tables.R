library(reshape2)

options(stringsAsFactors=F)

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
## Rm duplicate biopsies and keep samples with >=50 indels and >=30 SVs
analysis_sample_selection <- read.delim(paste0(base_dir,'/HMF_DR010_DR047/analysis/pancancer_overview/data/analysis_sample_selection.txt'))
selected_samples <- analysis_sample_selection[analysis_sample_selection$is_selected,'sample_id']

#--------- CHORD ---------#
chord_dir <- paste0(base_dir,'/HMF_DR010_DR047/scripts/chord_training/main/output/60_snvContext_svContext_allSigsRel_noSuppEvidence_customBlacklist3_noBoruta/seed015')
pred <- (function(){
   df <- read.delim(paste0(chord_dir,'/whole_hmf_dataset_probs.txt'))
   df <- cbind(sample=rownames(df),df)
   rownames(df) <- NULL
   #df <- df[df$sample %in% selected_samples,]
   return(df)
})()

hrd_samples <- pred$sample[pred$hrd>=0.5]


#--------- Determine MSI samples ---------#
contexts <- readRDS(paste0(base_dir,'/HMF_DR010_DR047/matrices/merged_contexts.rds'))

msi_samples <- (function(){
   contexts_indel_simple <- contexts$indel[,grep('[.]rep[.]',colnames(contexts$indel))]
   indel_load <- rowSums(contexts_indel_simple)
   names(indel_load)[indel_load>14000]
})()

#--------- Genotypes ---------#
l_m_diplotypes <- readRDS(paste0(base_dir,'/HMF_DR010_DR047/analysis/pancancer_overview/data/l_m_diplotypes.rds'))
rank_order_clust <- readRDS(paste0(base_dir,'/HMF_DR010_DR047/analysis/pancancer_overview/data/rank_order_clust.rds'))

## Used pre-subsetted and filtered diplotypes
diplotypes_hrd_hrGenes <- read.delim(paste0(base_dir,'/HMF_DR010_DR047/analysis/pancancer_overview/data/diplotypes_hrd_hrGenes.txt.gz'))
diplotypes_hrGenes <- read.delim(paste0(base_dir,'/HMF_DR010_DR047/analysis/pancancer_overview/data/diplotypes_hrGenes.txt.gz'))

diplotypes_raw <- (function(){
   path <- '/Users/lnguyen/Desktop/hmf_gene_diplotypes_max.txt.gz'
   
   if(file.exists(path)){
      df <- read.delim(path)
   } else {
      df <- read.delim(paste0(base_dir,'/HMF_DR010_DR047/scripts/annotate_genes/hmf_gene_diplotypes_max.txt.gz'))
      write.table(df,gzfile(path), sep='\t', quote=F, row.names=F)
   }
   
   return(df)
})()

diplotypes <- diplotypes_raw[diplotypes_raw$sample %in% selected_samples,]

#subset(diplotypes, !(sample %in% hrd_samples) & hit_score==10)

#--------- Other ---------#
metadata <- read.delim(paste0(base_dir,'/HMF_DR010_DR047/metadata/DR-047-update1_metadata.txt'))
hmf_sample_ids <- structure(metadata$hmf_id, names=metadata$sample_id)

ann_data_dir <- '/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD_data/HMF_DR010_DR047/vcf_subset/'

vsub <- function(v, subs){
   #v=pd_rel$L1
   #subs=c(hrd='HRD',hrp='HRP')
   
   subs[match(v, names(subs))]
}


####################################################################################################
# Full diplotypes                                                                                  #
####################################################################################################

#========= Main function =========#
getFullDiplotypes <- function(diplotypes, ann.data.dir=ann_data_dir){
   # diplotypes=diplotypes_hrd_hrGenes
   # ann.data.dir=ann_data_dir
   
   diplotypes_x_sample <- split(diplotypes, diplotypes$sample)
   
   getSampleFullDiplotypes <- function(sample_diplotypes, sample_name=NULL){
      # sample_diplotypes=diplotypes_x_sample[[1]]
      
      if(is.null(sample_name)){ sample_name <- sample_diplotypes$sample[1] }
      
      ## Read full diplotypes and subset by genes of interest (for speed)
      diplotypes_full <- read.delim(paste0(ann_data_dir,'/',sample_name,'/gene_statuses/gene_diplotypes.txt.gz'))
      diplotypes_full <- diplotypes_full[diplotypes_full$ensembl_gene_id %in% sample_diplotypes$ensembl_gene_id,]
      
      diplotypes_full_sel <- lapply(1:nrow(sample_diplotypes), function(i){
         
         #i=2
         r <- sample_diplotypes[i,]
         
         out <- diplotypes_full[
            diplotypes_full$ensembl_gene_id == r$ensembl_gene_id
            & diplotypes_full$diplotype_origin == r$diplotype_origin
            & diplotypes_full$hit_type == r$hit_type
            
            & diplotypes_full$a1 == r$a1
            & diplotypes_full$a2 == r$a2
            
            & diplotypes_full$a1.max_score == r$a1.max_score_old
            & diplotypes_full$a2.max_score == r$a2.max_score_old
            ,]
         
         ## Set pon variants to 0 / 'none'
         if(r$a1.in_pon){
            out <- within(out,{
               a1.chrom <- 'none'
               a1.pos <- 0
               a1.hgvs_c <- 'none'
               a1 <- 'none'
               a1.max_score <- 0
               a1.max_score_origin <- 'none'
            })
         }
         
         if(r$a2.in_pon){
            out <- within(out,{
               a2.chrom <- 'none'
               a2.pos <- 0
               a2.hgvs_c <- 'none'
               a2 <- 'none'
               a2.max_score <- 0
               a2.max_score_origin <- 'none'
            })
         }
         
         ## Deal with 'none', 'none' (a1, a2)
         out <- unique(out)
         
         ## Add sample name
         out <- cbind(sample=sample_name, out)
         
         return(out)
      })
      
      do.call(rbind, diplotypes_full_sel)
   }
   
   counter <- 0
   #pb <- txtProgressBar(max=length(diplotypes_x_sample), style=3)
   l_diplotypes_full_ss <- lapply(diplotypes_x_sample, function(i){
      counter <<- counter + 1
      sample_name <- names(diplotypes_x_sample)[counter]
      message('Processing [',counter,'/',length(diplotypes_x_sample),']: ', sample_name)
      #setTxtProgressBar(pb, value=counter)
      
      getSampleFullDiplotypes(i, sample_name = sample_name)
   })
   
   out <- do.call(rbind,l_diplotypes_full_ss)
   rownames(out) <- NULL
   return(out)
}


#========= Formatting functions =========#
formatFullDiplotyes <- function(diplotypes){
   #diplotypes=diplotypes_hrd_hrGenes_full
   #colnames(diplotypes)
   
   sel_cols <- c(
      sample=NA,
      ensembl_gene_id=NA,
      hgnc_symbol=NA,
      hit_score='biallelic_p_score',
      #hit_score_boosted,
      #diplotype_origin,
      hit_type='biallelic_hit_type',
      a1.chrom=NA,   
      a1.pos=NA,
      a1.hgvs_c=NA,
      a1='a1.event',
      a1.max_score='a1.p_score',
      a1.max_score_origin='a1.p_score_evidence',
      a2.chrom=NA,
      a2.pos=NA,
      a2.hgvs_c=NA,      
      a2='a2.event',
      a2.max_score='a2.p_score',
      a2.max_score_origin='a2.p_score_evidence'
   )
   
   sel_cols[is.na(sel_cols)] <- names(sel_cols)[is.na(sel_cols)]
   
   out <- diplotypes[names(sel_cols)]
   colnames(out) <- sel_cols
   
   return(out)
}

addSampleIndex <- function(diplotypes){
   cbind(
      sample_index=match(diplotypes$sample, names(rank_order_clust$clusters)),
      diplotypes
   )
}

addClusterNumber <- function(diplotypes){
   cbind(
      cluster=rank_order_clust$clusters[ match(diplotypes$sample, names(rank_order_clust$clusters)) ],
      diplotypes
   )
}

sampleNames2HmfIds <- function(v){ vsub(v, hmf_sample_ids) }

renameDeepDelTypes <- function(v){
   v[v=='full_gene_loss'] <- 'deep_deletion'
   v[v=='trunc'] <- 'deep_deletion'
   return(v)
}

getHrdProbs <- function(v){ 
   out <- pred[ match(v, pred$sample), c('hrd','BRCA1','BRCA2') ] 
   colnames(out) <- c('prob_hrd','prob_brca1_type_hrd','prob_brca2_type_hrd')
   return(out)
}

#========= HRD =========#
diplotypes_hrd_hrGenes_full_path <- paste0(base_dir,'/HMF_DR010_DR047/analysis/pancancer_overview/data/diplotypes_hrd_hrGenes_full.txt.gz')

if(!file.exists(diplotypes_hrd_hrGenes_full_path)){
   diplotypes_hrd_hrGenes_full <- getFullDiplotypes(diplotypes_hrd_hrGenes)
   write.table(
      diplotypes_hrd_hrGenes_full, gzfile(diplotypes_hrd_hrGenes_full_path),
      sep='\t',row.names=F,quote=F
   )
} else {
   diplotypes_hrd_hrGenes_full <- read.delim(diplotypes_hrd_hrGenes_full_path)
}

diplotypes_hrd_hrGenes_full_export <- (function(){
   
   df <- diplotypes_hrd_hrGenes_full
   
   ## Reorder based on rank order clustering
   df <- do.call(rbind,{
      split(df, df$sample)[
         names(rank_order_clust$clusters) 
      ]
   })
   
   df <- formatFullDiplotyes(df)
   
   df <- addClusterNumber(df)
   df <- addSampleIndex(df)
   
   df$a1.event <- renameDeepDelTypes(df$a1.event)
   df$a2.event <- renameDeepDelTypes(df$a2.event)
   df$biallelic_hit_type <- renameDeepDelTypes(df$biallelic_hit_type)
   
   df <- cbind(df,getHrdProbs(df$sample))
   df$sample <- sampleNames2HmfIds(df$sample)
   
   rownames(df) <- NULL
   
   df <- df[order(df$sample_index, df$hgnc_symbol),]
   
   return(df)
})()

#========= HRP samples with high hit_score =========#
diplotypes_hrp_hrGenes <- subset(diplotypes_hrGenes, !(sample %in% hrd_samples))

diplotypes_hrp_hrGenes_ss <- subset(diplotypes_hrp_hrGenes, hit_score>=9)

diplotypes_hrp_hrGenes_full_path <- paste0(base_dir,'/HMF_DR010_DR047/analysis/pancancer_overview/data/diplotypes_hrp_hrGenes_full.txt.gz')
if(!file.exists(diplotypes_hrp_hrGenes_full_path)){
   diplotypes_hrp_hrGenes_full <- getFullDiplotypes(diplotypes_hrp_hrGenes_ss)
   write.table(
      diplotypes_hrp_hrGenes_full, gzfile(diplotypes_hrp_hrGenes_full_path),
      sep='\t',row.names=F,quote=F
   )
} else {
   diplotypes_hrp_hrGenes_full <- read.delim(diplotypes_hrp_hrGenes_full_path)
}



diplotypes_hrp_hrGenes_full_export <- (function(){
   df <- diplotypes_hrp_hrGenes_full
   
   df <- formatFullDiplotyes(df)
   
   df$a1.event <- renameDeepDelTypes(df$a1.event)
   df$a2.event <- renameDeepDelTypes(df$a2.event)
   df$biallelic_hit_type <- renameDeepDelTypes(df$biallelic_hit_type)
   
   df <- cbind(df,getHrdProbs(df$sample))
   df$has_msi <- df$sample %in% msi_samples
   
   df$sample <- sampleNames2HmfIds(df$sample)
   
   return(df)
})()

####################################################################################################
# Get VUS's                                                                                        #
####################################################################################################

#========= Narrow down VUS's =========#
splitDiplotypeMatrix <- function(m){
   l <- list(
      a1=m[,grep('_1$',colnames(m))],
      a2=m[,grep('_2$',colnames(m))]
   )
}

meltDiplotypeMatrix <- function(l_m_diplotypes){
   meltDiplotypeMatrixAllele <- function(allele='a2'){
      counter <- 0
      l <- lapply(l_m_diplotypes, function(i){
         counter <<- counter + 1
         df <- as.data.frame(splitDiplotypeMatrix(i)[[allele]])
         df <- cbind(sample=rownames(df), df)
         df_melt <- melt(df, 'sample')
         
         df_melt$variable <- sapply(strsplit(as.character(df_melt$variable),'_'),`[`,1)
         colnames(df_melt)[2] <- 'hgnc_symbol'
         colnames(df_melt)[3] <- names(l_m_diplotypes)[counter]
         return(df_melt)
      })
      
      Reduce(function(x,y){ merge(x,y,sort=F) },l)
   }
   
   list(
      a1=meltDiplotypeMatrixAllele('a1'),
      a2=meltDiplotypeMatrixAllele('a2')
   )
}

getDiplotypesWithNovelVus <- function(l_m_diplotypes, rank_order_clust, diplotypes){
   
   #--------- Make filters ---------#
   ## Convert diplotype matrix back to diplotype table
   l_m_diplotypes_melt <- meltDiplotypeMatrix(l_m_diplotypes)
   
   df_selection <- data.frame(
      sample=l_m_diplotypes_melt$a1$sample,
      hgnc_symbol=l_m_diplotypes_melt$a1$hgnc_symbol,
      a1=l_m_diplotypes_melt$a1$eff,
      a2=l_m_diplotypes_melt$a2$eff,
      a1.max_score=l_m_diplotypes_melt$a1$score,
      a2.max_score=l_m_diplotypes_melt$a2$score,
      a2.max_score_origin=l_m_diplotypes_melt$a2$max_score_origin
   )
   
   df_selection$has_new_pathogenic <- with(df_selection,{
      # a1.max_score >= 5 & a2.max_score >= 3 & 
      #    (a2.max_score < 5 | (a2=='frameshift' & a2.max_score_origin!='known'))
      a1.max_score >= 5 & a2.max_score >= 3 & a2.max_score < 5
   })
   
   ## Get cluster per sample
   sample_clusters <- rank_order_clust$clusters
   cluster_names <- rank_order_clust$cluster_names
   
   for(i in cluster_names){
      sample_clusters[sample_clusters==i] <- names(cluster_names)[i]
   }
   
   df_selection$cluster <- sample_clusters
   
   ## Subset for potentially new pathogenic variants
   df_selection <- df_selection[df_selection$has_new_pathogenic,]
   
   
   #--------- Get diplotypes ---------#
   diplotypes_ss <- do.call(rbind,lapply(1:nrow(df_selection), function(i){
      row <- df_selection[i,]
      
      variant_type_selection <- row$a2
      
      if(variant_type_selection=='frameshift'){ variant_type_selection <- 'frameshift_variant' }
      else if(variant_type_selection=='essential_splice_variant'){ variant_type_selection <- c('splice_acceptor_variant','splice_donor_variant') }
      else if(variant_type_selection=='nonsense'){ variant_type_selection <- c('start_lost','stop_gained','stop_lost') }
      else if(variant_type_selection=='missense'){ variant_type_selection <- 'missense_variant' }
      
      diplotypes[
         diplotypes$sample==row$sample 
         & diplotypes$hgnc_symbol==row$hgnc_symbol
         & diplotypes$a2 %in% variant_type_selection
         & !(diplotypes$diplotype_origin %in% c('germ_som','som_som'))
         ,]
   }))
   
   diplotypes_ss$cluster <- sample_clusters[ match(diplotypes_ss$sample,names(sample_clusters)) ]
   diplotypes_ss <- diplotypes_ss[diplotypes_ss$cluster==diplotypes_ss$hgnc_symbol,]
   
   return(diplotypes_ss)
}

diplotypes_vus <- getDiplotypesWithNovelVus(l_m_diplotypes, rank_order_clust, diplotypes_hrd_hrGenes)


#========= Get variants from raw data =========#
getVusFromDiplotype <- function(
   diplotypes,
   ann.data.dir='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD_data/HMF_DR010_DR047/vcf_subset/'
){
   # diplotypes=diplotypes_vus
   # ann.data.dir='/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD_data/HMF_DR010_DR047/vcf_subset/'
   
   #--------- Search through raw data ---------#
   l <- with(diplotypes,{
      Map(function(sample, ensembl_gene_id, diplotype_origin, a2.max_score){
         
         a2_origin <- strsplit(diplotype_origin,'_')[[1]][2]
         ann_data_path <- paste0(ann.data.dir,'/',sample,'/gene_statuses/mut_profile_',a2_origin,'.txt.gz')
         
         print(ann_data_path)
         
         df <- read.delim(ann_data_path)
         df_ss <- df[df$ensembl_gene_id==ensembl_gene_id & df$max_score==a2.max_score,]
         df_ss$mut_origin <- a2_origin
         
         return(df_ss)
         
      }, sample, ensembl_gene_id, diplotype_origin, a2.max_score)
   })
   
   #--------- Merge and format columns/values ---------#
   out <- do.call(rbind, l)
   sel_cols <- c(
      'chrom'='chrom',
      'pos'='pos',
      'ref'='ref',
      'alt'='alt',
      'hgvs_c'='hgvs_c',
      'snpeff_eff'='variant_type',
      'mut_origin'='mutation_origin',
      'hgnc_symbol'='hgnc_symbol',
      'ensembl_gene_id'='ensembl_gene_id',
      'clinvar_sig'='clinvar_annotation'
   )
   out <- out[names(sel_cols)]
   
   out <- within(out,{
      clinvar_sig[is.na(clinvar_sig)] <- ''
      mut_origin[mut_origin=='germ'] <- 'germline'
      mut_origin[mut_origin=='som'] <- 'somatic'
   })
   
   colnames(out) <- sel_cols
   
   out <- cbind(
      sample=sapply(strsplit(rownames(out),'[.]'),`[[`,1),
      out
   )
   rownames(out) <- NULL
   return(out)
}

mut_profile_vus <- getVusFromDiplotype(diplotypes_vus)

mut_profile_vus_export <- (function(){
   df <- mut_profile_vus
   df <- addClusterNumber(df)
   df <- addSampleIndex(df)
   df$sample <- sampleNames2HmfIds(df$sample)
   
   df <- df[order(df$sample_index),]
   
   return(df)
})()

write.table(
   mut_profile_vus, paste0(base_dir,'/HMF_DR010_DR047/analysis/pancancer_overview/supp_tables/hmf_hrd_novel_pathogenic.txt'),
   sep='\t',row.names=F,quote =F
)

####################################################################################################
# Export xlsx                                                                                      #
####################################################################################################

write.tsv <- function(...){
   write.table(..., sep='\t',row.names=F,quote=F)
}

write.tsv(
   diplotypes_hrd_hrGenes_full_export,
   paste0(base_dir,'/HMF_DR010_DR047/analysis/pancancer_overview/supp_tables/diplotypes_hrd_samples.txt')
)

write.tsv(
   mut_profile_vus_export,
   paste0(base_dir,'/HMF_DR010_DR047/analysis/pancancer_overview/supp_tables/likely_pathogenic_vus.txt')
)

write.tsv(
   diplotypes_hrp_hrGenes_full_export,
   paste0(base_dir,'/HMF_DR010_DR047/analysis/pancancer_overview/supp_tables/diplotypes_hrp_samples_hi_hit_score.txt')
)








