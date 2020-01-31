library(ggplot2)
library(cowplot); theme_set(theme_grey())

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

main <- function(output_dir){
  annotation <- read.delim(paste0(output_dir,'/annotation.txt'))
  
  #========= Extract CV predictions =========#
  gatherCvPred <- function(dir){

    fold_dirs <- list.dirs(paste0(dir,'/cv/'), recursive=F)
    l <- lapply(paste0(fold_dirs,'/rf_out.rds'), function(i){
      readRDS(i)$pred
    })
    
    df <- do.call(rbind, l)
    df <- cbind(sample=rownames(df), df)
    rownames(df) <- NULL
    return(df)
  }
  
  cv_pred_path <- paste0(output_dir,'/cv_pred.rds')
  
  if(file.exists(cv_pred_path)){
    l_pred <- readRDS(cv_pred_path)
  } else {
    dirs <- dir(paste0(output_dir,'/replicates/'), recursive=F, full.names=T, pattern='seed_')
    
    counter <- 0
    l_pred <- lapply(dirs, function(i){
      counter <<- counter+1
      message('Processing [',counter,']')
      
      df <- gatherCvPred(i)
      df <- cbind(set_id=counter, set_name=basename(i),df)
      df <- df[order(df$sample),]
      df$hrd <- df$BRCA1 + df$BRCA2
      df$is_hrd <- df$hrd >= 0.5
      
      return(df)
    })
    names(l_pred) <- basename(dirs)
    saveRDS(l_pred, cv_pred_path)
  }
  
  
  #========= Counts times predicted as HRD =========#
  cv_pred_stats <- (function(){
    gatherCols <- function(col){
      m <- do.call(cbind, lapply(l_pred, function(i){ i[,col] }))
      rownames(m) <- l_pred[[1]]$sample
      return(m)
    }
    
    n_is_hrd <- rowSums(gatherCols('is_hrd'))
    
    hrd <- gatherCols('hrd')
    prob_mean <- rowMeans(hrd)
    prob_sd <- apply(hrd,1,sd)
    
    out <- data.frame(n_is_hrd,prob_mean,prob_sd)
    out <- cbind(sample=rownames(out),out)
    
    #c('sample','response','in_training_set')
    annotation_ss <- annotation[match(out$sample, annotation$sample),]
    out <- merge(out, annotation_ss, by='sample')
    out <- out[order(out$prob_mean,decreasing=T),]
    rownames(out) <- NULL
    return(out)
  })()
  
  #========= Plot =========#
  MIN_PROP_IS_HRD_BRCA <- 0.6
  MAX_PROP_IS_HRD_NONE <- 0.4
  
  plotCvPredStats <- function(
    df, export.path=NULL, pdf.dim=NULL, 
    min.prop.is.hrd.brca=MIN_PROP_IS_HRD_BRCA,
    max.prop.is.hrd.none=MAX_PROP_IS_HRD_NONE
  ){
    #df <- subset(cv_pred_stats, response!='none' & in_training_set)
    #df <- cv_pred_stats
    
    df$sample <- factor(df$sample, df$sample)
    
    if(is.null(pdf.dim)){
      pdf.dim <- c(nrow(df)/10, 8)
    }
    
    ## plot parameters
    color_pal <- c(BRCA1='red',BRCA2='blue',none='black')
    
    plot_theme <- theme_bw() + theme(
      #plot.title=element_text(hjust=0.5),
      axis.text.x=element_text(angle=90,vjust=0.5, hjust=1),
      axis.title.x=element_blank(),
      panel.border=element_rect(fill=NA)
    )
    
    ## main
    plots <- list()
    
    upper_hline <- max(df$n_is_hrd) * MIN_PROP_IS_HRD_BRCA
    lower_hline <- max(df$n_is_hrd) * MAX_PROP_IS_HRD_NONE
    plots$n_is_hrd <- ggplot(df, aes(sample, n_is_hrd, color=response)) + 
      geom_point() + 
      
      geom_hline(yintercept=upper_hline, linetype='dashed') +
      annotate(
        "text", 
        x=1, hjust=0, y=upper_hline, vjust=-1,
        label=sprintf("Min times HRD = %s: 'BRCA1','BRCA2' class whitelist", upper_hline)
      ) +
      
      geom_hline(yintercept=lower_hline, linetype='dashed') +
      annotate(
        "text", 
        x=1, hjust=0, y=lower_hline, vjust=1.5,
        label=sprintf("Max times HRD = %s: 'none' class whitelist", lower_hline)
      ) +
      
      ggtitle(paste0('N = ',nrow(df))) +
      ylab('n times HRD') +
      scale_color_manual(values=color_pal) +
      plot_theme
    
    plots$prob_mean <- ggplot(df, aes(sample, prob_mean, color=response)) + 
      geom_point() + ylab('Mean HRD prob.') +
      geom_errorbar(aes(ymin=prob_mean-prob_sd, ymax=prob_mean+prob_sd)) +
      scale_color_manual(values=color_pal) +
      plot_theme
    
    out <- plot_grid(plotlist=plots, align='v', ncol=1, rel_heights=c(1.1,1))
    if(is.null(export.path)){
      return(out)
    } else {
      pdf(export.path,pdf.dim[1],pdf.dim[2])
      plot(out)
      dev.off()
    }
  }
  
  #plotCvPredStats(subset(cv_pred_stats, response!='none' & in_training_set))
  
  plotCvPredStats(
    subset(cv_pred_stats, response!='none' & in_training_set),
    export.path=paste0(output_dir,'/cv_pred_stats.pdf')
  )
  
  plotCvPredStats(
    cv_pred_stats,
    export.path=paste0(output_dir,'/cv_pred_stats_full.pdf')
  )
  
  
  #========= Make blacklist =========#
  custom_blacklist <- list(
    brca=subset(cv_pred_stats, response %in% c('BRCA1','BRCA2') & n_is_hrd <= max(n_is_hrd)* MIN_PROP_IS_HRD_BRCA)$sample,
    none=subset(cv_pred_stats, response == 'none' & n_is_hrd >= max(n_is_hrd) * MAX_PROP_IS_HRD_NONE)$sample
  )
  
  saveRDS(custom_blacklist, paste0(output_dir,'/custom_blacklist.rds'))
}

fold_dirs <- list.dirs(paste0(base_dir,'/HMF_DR010_DR047/scripts/chord_training/main/output/60_CV_snvContext_svContext_allSigsRel_noSuppEvidence_noBoruta/'), recursive=F)
for(i in fold_dirs){
  message(paste0('Processing: ', basename(i)))
  main(i)
}



