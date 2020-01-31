library(mutSigExtractor)
library(ggplot2)
library(scales)
library(grid)
library(ggpubr)

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

#========= Misc functions =========#
forceDfOrder <- function(df){
   as.data.frame(lapply(df, function(i){
      if(!is.numeric(i)){ i <- factor(i, unique(i)) }
      return(i)
   }))
}

#========= Load data =========#
analysis_sample_selection <- read.delim(paste0(base_dir,'/HMF_DR010_DR047/analysis/pancancer_overview/data/analysis_sample_selection.txt'))
selected_samples <- analysis_sample_selection[analysis_sample_selection$is_selected,'sample_id']

#--------- Contexts ---------#
contexts <- list()

contexts$abs <- readRDS(paste0(base_dir,'/HMF_DR010_DR047/matrices/merged_contexts.rds'))
contexts$abs <- transformContexts(contexts$abs, simplify.types=c('snv','indel'), export.list=T)

contexts$rel <- lapply(contexts$abs, function(i){ 
   t(apply(i,1,function(j){ j/sum(j) }) )
})

#--------- Annotation ---------#
pred_hmf <- read.delim(paste0(base_dir,'/HMF_DR010_DR047/analysis/chord_perf/pred_tsv/pred_hmf.txt'))

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

#sel_genes <- c('BRCA2','BRCA1','RAD51C','PALB2')
sel_genes <- c('BRCA2','BRCA1')
diplotypes_ss <- diplotypes[diplotypes$hgnc_symbol %in% sel_genes,]

diplotypes_split <- split(diplotypes_ss, diplotypes_ss$sample)

gene_defs <- unlist(lapply(diplotypes_split, function(i){
   gene_def_init <- 'none'
   gene_def <- i$hgnc_symbol[i$hit_score==10 & !(i$diplotype_origin %in% c('germ_som','som_som'))]
   if(length(gene_def)!=0){
      gene_def_init <- gene_def[1]
   }
   return(gene_def_init)
}))

annotation <- data.frame(
   sample=names(gene_defs),
   response=gene_defs,
   row.names=NULL
)

#========= Plot  =========#
#--------- Prep data ---------#
df <- pred_hmf[c('sample','p_chord')]
df$p_rank <- order(df$p_chord)
df <- merge(df,annotation,sort=F)

df$response2 <- df$response
#df$response2[df$custom_blacklist_checked] <- 'none'

df$del.mh <- contexts$abs$indel[match(df$sample, rownames(contexts$abs$indel)),'del.mh'] 
df$del.mh_rel <- contexts$rel$indel[match(df$sample, rownames(contexts$rel$indel)),'del.mh'] 

df <- cbind(
   df,
   (function(){
      m <- as.data.frame(contexts$abs$indel[,c('del.rep','ins.rep')] )
      m$indel.rep <- apply(m,1,sum)
      m$indel.all <- rowSums(contexts$abs$indel)
      m$indel.rep_rel <- m$indel.rep / m$indel.all
      
      m[match(df$sample, rownames(m)),]
   })()
)

df <- df[order(df$response, decreasing=T),] ## Force plotting of BRCA1/2 deficient samples last
df <- as.data.frame(lapply(df, function(i){
   if(!is.numeric(i)){ i <- factor(i, unique(i)) }
   return(i)
}))

## Plot points with gene deficiencies last
df <- do.call(rbind,split(df,df$response2!='none'))

df$hr_status <- ifelse(df$p_chord >= 0.5, 'hrd', 'hrp')
df$has_msi <- df$indel.rep>=14000
df$is_fn_hrp <- with(df,{ hr_status=='hrp' & response!='none' })

df_ss <- subset(df, is_fn_hrp)

## Ordering
df <- forceDfOrder(df)

gene_colors <- c(BRCA1='#f58225', BRCA2='#69439d', RAD51C='green', PALB2='magenta', none='grey')
df$response2 <- factor(df$response2, names(gene_colors))

#subset(df, indel.rep>14000 & p_chord<0.5 & response!='none')$sample


#--------- indel.rep vs. HRD score ---------#
plot_p_chord_vs_msi <- ggplot(df, aes(indel.rep, p_chord, color=response2)) + 
   geom_point() +
   geom_vline(xintercept=14000, linetype='dotted') +
   annotate('text',x=14000*1.1, y=0.5, label='MSI positive:\nindel.rep > 14000', hjust=0) +
   
   scale_color_manual(
      values=gene_colors,
      breaks=names(gene_colors)
   ) +
   scale_x_log10(labels=comma) +
   labs(y='Probability of HRD', x='indel.rep (abs. contrib.)', color='Gene deficiency') +
   #guides(guide_legend(reverse=F)) +
   theme_bw()
   


pdf(paste0(base_dir,'/HMF_DR010_DR047/analysis/chord_perf/chord_main/p_chord_vs_indel_rep.pdf'), 8, 5)
grid.draw(plot_p_chord_vs_msi)
dev.off()






