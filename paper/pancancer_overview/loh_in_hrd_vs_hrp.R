library(ggplot2)
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
#--------- Other data ---------#
chord_dir <- paste0(base_dir,'/HMF_DR010_DR047/scripts/chord_training/main/output/60_snvContext_svContext_allSigsRel_noSuppEvidence_customBlacklist3_noBoruta/seed015/')

#--------- CHORD predictions ---------#
## Rm duplicate biopsies and keep samples with >=50 indels
analysis_sample_selection <- read.delim(paste0(base_dir,'/HMF_DR010_DR047/analysis/pancancer_overview/data/analysis_sample_selection.txt'))
selected_samples <- analysis_sample_selection[analysis_sample_selection$is_selected,'sample_id']

pred <- (function(){
   df <- read.delim(paste0(chord_dir,'/whole_hmf_dataset_probs.txt'))
   df <- cbind(sample=rownames(df),df)
   rownames(df) <- NULL
   df <- df[df$sample %in% selected_samples,]
   return(df)
})()

hrd_samples <- pred$sample[pred$hrd >= 0.5]
hrp_samples <- pred$sample[pred$hrd < 0.5]

n_hrd <- length(hrd_samples)
n_hrp <- length(hrp_samples)

#--------- Diplotypes ---------#
diplotypes_filt <- read.delim(paste0(base_dir,'/HMF_DR010_DR047/analysis/pancancer_overview/data/diplotypes_filt.txt.gz'))

diplotypes_filt_ss <- diplotypes_filt[
   diplotypes_filt$hgnc_symbol %in% c('BRCA1','BRCA2','RAD51C','PALB2')
,]

#========= Find presence of LOH in any of the 4 genes =========#
calcCnvFreq <- function(alteration, samples){
   length(
      with(diplotypes_filt_ss,{
         unique( sample[sample %in% samples & a1 %in% alteration] )
      })
   )
}

freqs <- list(
   hrd_loh=calcCnvFreq('loh',hrd_samples),
   hrp_loh=calcCnvFreq('loh',hrp_samples),
   
   hrd_dd=calcCnvFreq(c('full_gene_loss','trunc'),hrd_samples),
   hrp_dd=calcCnvFreq(c('full_gene_loss','trunc'),hrp_samples)
)

m_conting <- matrix(
   c( freqs$hrd_loh, n_hrd, freqs$hrp_loh, n_hrp ),
   nrow=2, byrow=T,
   dimnames=list(c('hrd','hrp'),c('has_loh','total'))
)
fisher <- fisher.test(m_conting, alternative='greater')

#========= Plotting =========#
pd <- data.frame(
   has_loh=c(freqs$hrd_loh, freqs$hrp_loh),
   not_has_loh=c(n_hrd-freqs$hrd_loh, n_hrp-freqs$hrp_loh)
)
rownames(pd) <- c('HRD','HRP')

pd_prop <- as.data.frame(pd/c(n_hrd,n_hrp))
pd_prop$hr_status <- rownames(pd_prop)

pd_prop$labels <- c(
   paste0(freqs$hrd_loh,' / ',n_hrd),
   paste0(freqs$hrp_loh,' / ',n_hrp)
)

pd_prop$label_pos <- pd_prop$has_loh/2

p <- ggplot(pd_prop, aes(hr_status, has_loh)) + 
   geom_bar(stat='identity', fill=c('#f09aa3','#b0e1a1'), color='black', size=0.3) +
   geom_text(mapping=aes(y=label_pos, label=labels), vjust=0.5) +
   
   geom_segment(aes(x=0.9, xend=2.1, y=0.95, yend=0.95), size=0.25) +
   annotate('text', x=1.5, y=1.05, label=paste0('p=',signif(fisher$p.value,4),"\nFisher's exact test")) +
   
   scale_y_continuous(labels=function(x){ paste0(x*100,'%') },limits=c(0, 1.1)) +
   
   labs(
      y='% with LOH in\nBRCA1/2, RAD51C or PALB2',
      x='HR status (CHORD prediction)'
   ) + 
   
   theme_bw() +
   theme(
      axis.text=element_text(size=10.5)
   )

pdf(paste0(base_dir,'/HMF_DR010_DR047/analysis/pancancer_overview/plots/loh_in_hrd_vs_hrp.pdf'), 4,4)
plot(p)
dev.off()





