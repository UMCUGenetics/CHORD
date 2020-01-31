library(mutSigExtractor)
library(ggpubr)
#library(reshape2)
#library(grid)
#library(ggplot2)
#library(cowplot)#; theme_set(theme_grey())

#library(biomaRt)

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

#--------- CHORD predictions ---------#
pred <- (function(){
   df <- read.delim(paste0(base_dir,'/HMF_DR010_DR047/scripts/chord_training/main/output/60_snvContext_svContext_allSigsRel_noSuppEvidence_customBlacklist3_noBoruta/seed015/whole_hmf_dataset_probs.txt'))
   df <- cbind(sample=rownames(df),df)
   rownames(df) <- NULL
   return(df)
})()

hrd_samples <- pred[pred$hrd>=0.5,'sample']

#--------- Read and cache diplotype table ---------#
## All cancer genes
diplotypes <- (function(){
   path <- '/Users/lnguyen/Desktop/hmf_gene_diplotypes_max.txt.gz'
   
   if(file.exists(path)){
      df <- read.delim(path)
   } else {
      df <- read.delim(paste0(base_dir,'/HMF_DR010_DR047/scripts/annotate_genes/hmf_gene_diplotypes_max.txt.gz'))
      write.table(df,gzfile(path), sep='\t', quote=F, row.names=F)
   }
   
   return(df)
})()

sel_genes <- c('BRCA1','BRCA2','RAD51C','PALB2')
diplotypes_2fs <- diplotypes[
   with(diplotypes,{ hgnc_symbol %in% sel_genes & n_frameshifts>=2 & hit_score==10 & a1=='loh' })
,]
#diplotypes_2fs[,c('sample','ensembl_gene_id','hgnc_symbol','a1','a2','n_frameshifts','fs_lengths','fs_origins','fs_positions','reading_frame')]

#subset(diplotypes_2fs, sample %in% selected_samples)

diplotypes_ss <- diplotypes[
   diplotypes$hgnc_symbol %in% sel_genes 
   & diplotypes$sample %in% selected_samples
   ,]

diplotypes_split <- split(diplotypes_ss, diplotypes_ss$sample)

gene_defs <- unlist(lapply(diplotypes_split, function(i){
   gene_def_init <- 'none'
   gene_def <- i$hgnc_symbol[i$hit_score==10 & !(i$diplotype_origin %in% c('germ_som','som_som'))]
   if(length(gene_def)!=0){
      gene_def_init <- gene_def[1]
   }
   return(gene_def_init)
}))

#--------- Clonal/subclonal indels ---------#
contexts_split_path <- paste0(base_dir,'/HMF_DR010_DR047/matrices_subclonal/merged_indel_contexts.rds')

if(file.exists(contexts_split_path)){
   contexts_split <- readRDS(contexts_split_path)
} else {
   contexts_split <- (function(){
      paths <- list.files(paste0(base_dir,'/HMF_DR010_DR047/matrices_subclonal/indel/'), full.names=T)
      contexts <- lapply(paths, read.delim)
      names(contexts) <- basename(paths)
      
      l <- list(
         clonal=contexts[grep('_CLONAL.txt$',names(contexts))],
         subclonal=contexts[grep('_SUBCLONAL.txt$',names(contexts))]
      )
      
      l <- lapply(l, function(i){
         m <- t(do.call(cbind,i))
         transformContexts(indel=m, simplify.types='indel')
      })
      
      return(l)
   })()
   
   saveRDS(contexts_split, paste0(base_dir,'/HMF_DR010_DR047/matrices_subclonal/merged_indel_contexts.rds'))
}


#========= Compare clonal/subclonal indels =========#
#contexts_split$subclonal[rownames(contexts_split$subclonal) %in% diplotypes_2fs$sample,]

## Remove samples with 100 or fewer subclonal indels
contexts_split_filt_rel <- lapply(contexts_split, function(i){
   m <- i[rowSums(contexts_split$subclonal) >= 100,]
   m/rowSums(m)
})

## Select del mh
del_mh_counts <- data.frame(
   clonal=contexts_split_filt_rel$clonal[,'del.mh'], 
   subclonal=contexts_split_filt_rel$subclonal[,'del.mh']
)

del_mh_counts <- cbind(sample=rownames(del_mh_counts),del_mh_counts)
rownames(del_mh_counts) <- NULL

## Add annotations
del_mh_counts$has_2fs <- del_mh_counts$sample %in% diplotypes_2fs$sample
del_mh_counts$fs_status <- ifelse(del_mh_counts$sample %in% diplotypes_2fs$sample, 'has_2fs', 'no_2fs')

del_mh_counts$is_hrd <- del_mh_counts$sample %in% hrd_samples
del_mh_counts$hr_status <- ifelse(del_mh_counts$sample %in% hrd_samples,'HRD','HRP')

del_mh_counts$group <- paste0(del_mh_counts$fs_status,'.',del_mh_counts$hr_status)

del_mh_counts <- del_mh_counts[del_mh_counts$sample %in% selected_samples,]

del_mh_counts$hrd_score <- pred$hrd[ match(del_mh_counts$sample, pred$sample) ]
del_mh_counts$response <- gene_defs[ match(del_mh_counts$sample, names(gene_defs)) ]

gene_colors <- c(BRCA1='#f58225', BRCA2='#69439d', RAD51C='magenta', PALB2='green', none='darkgrey')
del_mh_counts$response <- factor(del_mh_counts$response, names(gene_colors))

# subset(del_mh_counts, hr_status=='HRP' & hrd_score>0.1)
# ggplot(del_mh_counts, aes(subclonal, hrd_score)) + geom_point()

set.seed(1)
plot_subclonal_del_mh <- ggplot(del_mh_counts, aes(x=group, y=subclonal, group=group)) +
   facet_grid(~hr_status, scales='free_x', space='free_x') +
   geom_boxplot(outlier.colour=NA) +
   geom_jitter(aes(color=response), width=0.3) +
   labs(y='Proportion of del.mh\nin subclonal fraction') +
   scale_x_discrete(
      labels=c(
         has_2fs.HRD='Has secondary out-of-frame\nframeshift in HR gene', 
         no_2fs.HRD='Other HRD samples\n', 
         no_2fs.HRP='HRP'
      )
   ) +
   #scale_fill_distiller(palette='Spectral', name='HRD prob.') +
   scale_color_manual(values=gene_colors, breaks=sel_genes, name='Gene deficiency') +
   #scale_fill_manual(values=gene_colors, breaks=sel_genes, name='Gene deficiency') +
   theme_bw() +
   theme(
      plot.title=element_text(hjust=0.5),
      axis.title.x=element_blank()#,
      #legend.position='none'
   )


pdf(paste0(base_dir,'/HMF_DR010_DR047/analysis/pancancer_overview/plots/subclonal_del_mh.pdf'), 8, 3)
grid.draw(plot_subclonal_del_mh)
dev.off()


del_mh_counts_exp <- (function(){
   df <- del_mh_counts[order(del_mh_counts$has_2fs, del_mh_counts$hrd_score, decreasing=T),]
   sel_cols <- c(
      sample=NA,
      clonal='del.mh_clonal',
      subclonal='del.mh_subclonal',
      has_2fs='has_secondary_frameshift',
      is_hrd=NA,
      hrd_score='prob_hrd',
      response='gene_deficiency'
   )
   
   df <- df[,names(sel_cols)]
   colnames(df)[!is.na(sel_cols)] <- sel_cols[!is.na(sel_cols)]
   return(df)
})()

write.table(
   del_mh_counts_exp,
   paste0(base_dir,'/HMF_DR010_DR047/analysis/pancancer_overview/supp_tables/subclonal_del_mh.txt'),
   sep='\t',row.names=F,quote=F
)

