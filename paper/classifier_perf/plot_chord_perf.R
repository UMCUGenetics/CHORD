library(ggplot2)
library(ggpubr)
library(mutSigExtractor)
library(randomForest)
library(cowplot)
library(grid)
library(gridExtra)
library(mltoolkit)
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

#========= Misc functions =========#
forceDfOrder <- function(df){
  as.data.frame(lapply(df, function(i){
    if(!is.numeric(i)){ i <- factor(i, unique(i)) }
    return(i)
  }))
}


####################################################################################################
# Load data                                                                                        #
####################################################################################################

analysis_sample_selection <- read.delim(paste0(base_dir,'/HMF_DR010_DR047/analysis/pancancer_overview/data/analysis_sample_selection.txt'))
selected_samples <- analysis_sample_selection[analysis_sample_selection$is_selected,'sample_id']

#========= CHORD =========#
#--------- Annotation ---------#
annotation <- read.delim(paste0(base_dir,'/HMF_DR010_DR047/scripts/det_brca_status/hmf_brca_response.txt'))

annotation2 <- within(annotation,{
  response[ response != 'none' & has_msi ] <- 'none' ## Set msi samples to 'none'
  response[ response != 'none' & !(a1 %in% c('full_gene_loss','trunc','loh')) ] <- 'none' ## Set germ_som samples to 'none'
})

#--------- Final model ---------#
chord_dir <- paste0(base_dir,'/HMF_DR010_DR047/scripts/chord_training/main/output/60_snvContext_svContext_allSigsRel_noSuppEvidence_customBlacklist3_noBoruta/seed015')
chord <- readRDS(paste0(chord_dir,'/final/rf_out.rds'))$RF

#--------- 100x CV training set filter ---------#
cv100_pred <- readRDS(paste0(base_dir,'/HMF_DR010_DR047/scripts/chord_training/main/output/60_snvContext_svContext_allSigsRel_noSuppEvidence_customBlacklist3_noBoruta/cv100/cv_pred.rds'))

summarizeCv100Pred <- function(cv100_pred){
  gatherCols <- function(col){
    m <- do.call(cbind, lapply(cv100_pred, function(i){ i[,col] }))
    rownames(m) <- cv100_pred[[1]]$sample
    return(m)
  }
  
  n_is_hrd <- rowSums(gatherCols('is_hrd'))
  
  ## Calculate mean and sd for all prediction classes
  l <- lapply(c('hrd','BRCA1','BRCA2'), function(i){
    #i='hrd'
    probs <- gatherCols(i)
    prob_mean <- rowMeans(probs)
    prob_sd <- apply(probs,1,sd)
    
    df <- data.frame(mean=prob_mean, sd=prob_sd)
    colnames(df) <- paste0(i,'.',colnames(df))
    return(df)
  })
  
  out <- do.call(cbind, l)
  out <- cbind(sample=rownames(out), out); rownames(out) <- NULL
  out$n_is_hrd <- n_is_hrd[out$sample]
  
  annotation_ss <- annotation[match(out$sample, annotation$sample),]
  out <- merge(out, annotation_ss, by='sample')
  out <- out[order(out$hrd.mean),]
  out <- cbind(index=1:nrow(out),out)
  rownames(out) <- NULL
  return(out)
}

cv100_pred_stats <- summarizeCv100Pred(cv100_pred)

#--------- Outer CV ---------#
gatherOuterCv <- function(dir){
  fold_dirs <- list.dirs(dir, recursive=F)
  lapply(paste0(fold_dirs,'/final/rf_out.rds'), function(i){ readRDS(i) })
}

outer_cv <- gatherOuterCv(paste0(base_dir,'/HMF_DR010_DR047/scripts/chord_training/main/output/60_CV_snvContext_svContext_allSigsRel_noSuppEvidence_noBoruta/'))

gatherOuterCvPred <- function(outer_cv){
  
  l <- lapply(outer_cv, `[[`, 'pred')
  
  df <- as.data.frame(do.call(rbind, l))
  df$p_chord <- df$BRCA1 + df$BRCA2
  df <- df[order(df$response),]
  df <- cbind(sample=rownames(df), df); rownames(df) <- NULL
  return(df)
}

cv_pred_hmf <- gatherOuterCvPred(outer_cv)

#--------- HMF predictions ---------#
## Features
contexts_hmf <- readRDS(paste0(base_dir,'/HMF_DR010_DR047/matrices/merged_contexts.rds'))

## Main
chordPredictHmf <- function(model, contexts, annotation, mode='contexts'){
  if(mode=='contexts'){
    m <- transformContexts(contexts, simplify.types=c('snv','indel'), rel.types=c('snv','indel','sv'))
  } else if (mode=='contexts_abs') {
    m <- transformContexts(contexts, simplify.types=c('snv','indel'))
  } else if (mode=='sigs') {
    contexts_ss <- contexts
    contexts_ss$sv <- contexts_ss$sv[,!grepl('_0e00_',colnames(contexts_ss$sv))]
    m <- transformContexts(contexts_ss, simplify.types='indel', lsqnonneg.types=c('snv','sv'), rel.types=c('snv','indel','sv'))
  }
  
  df <- as.data.frame(predict(model, m, 'prob'))
  df$p_chord <- df$BRCA1 + df$BRCA2
  df <- cbind(sample=rownames(df), df)
  rownames(df) <- NULL
  
  sigs_rel2 <- cbind(sample=rownames(m),as.data.frame(m)); rownames(sigs_rel2) <- NULL
  
  df <- merge(df, sigs_rel2, by='sample')
  df$response <- annotation$response[match(df$sample, annotation$sample)]
  
  return(df)
}

pred_hmf <- chordPredictHmf(chord, contexts_hmf, annotation2)

detHrdType <- function(pred_hmf){
  which_max <- apply(pred_hmf[,c('BRCA1','BRCA2')],1,which.max)
  v <- c('BRCA1_type','BRCA2_type')[ which_max ]
  v[ pred_hmf$p_chord  < 0.5 ] <- 'none'
  return(v)
}

#--------- ICGC predictions ---------#
## Annotations
hrdetect_annotations <- (function(){
  paths <- list.files(paste0(base_dir,'/ICGC/analysis/hrdetect_output/'),pattern='hrdetect.*.txt',full.names=T)
  names(paths) <- sapply(strsplit(basename(paths),'[.]'),`[`,2)
  lapply(paths,read.delim)
})()

## Features
l_contexts_icgc <- (function(){
  dirs <- list.dirs(paste0(base_dir,'/ICGC/matrices/'),recursive=F)
  paths <- paste0(dirs,'/merged_contexts.rds')
  names(paths) <- gsub('-','_',basename(dirs))
  
  lapply(paths, function(i){ 
    readRDS(i)
  })
})()

l_contexts_icgc$BRCA_EU.HMF_pipeline <- readRDS(paste0(base_dir,'/BRCA_EU/matrices/BRCA_EU_hmfPipeline_2/merged_contexts.rds'))

## Main
chordPredictIcgc <- function(model, l_contexts_icgc, hrdetect_annotations, mode='contexts'){
  #--------- Model ---------#
  l_sigs_rel <- lapply(l_contexts_icgc, function(i){
    if(mode=='contexts'){
      m <- transformContexts(i, simplify.types=c('snv','indel'), rel.types=c('snv','indel','sv'))
    } else if (mode=='contexts_abs') {
      m <- transformContexts(i, simplify.types=c('snv','indel'))
    } else if (mode=='sigs') {
      contexts_ss <- i
      contexts_ss$sv <- contexts_ss$sv[,!grepl('_0e00_',colnames(contexts_ss$sv))]
      m <- transformContexts(contexts_ss, simplify.types='indel', lsqnonneg.types=c('snv','sv'), rel.types=c('snv','indel','sv'))
    }
    
    return(m)
  })
  
  #--------- Prediction ---------#
  l_pred <- lapply(l_sigs_rel, function(i){ 
    df <- as.data.frame(predict(model, i, 'prob'))
    df <- cbind(sample=rownames(df), df); rownames(df) <- NULL
    df$p_chord <- df$BRCA1 + df$BRCA2
    return(df)
  })
  
  #--------- Merge CHORD and HRDetect data ---------#
  l_sigs_rel2 <- lapply(l_sigs_rel, function(i){
    i <- as.data.frame(i)
    i <- cbind(sample=rownames(i),i)
    rownames(i) <- NULL
    return(i)
  })
  
  merged_data <- lapply(names(l_pred), function(i){
    pred=l_pred[[i]]
    
    cohort_name <- sapply(strsplit(i,'[.]'),`[[`,1)
    hrdetect.annotations=hrdetect_annotations[[cohort_name]]
    
    sigs.rel=cbind(sample=rownames(l_sigs_rel2[[i]]), as.data.frame(l_sigs_rel2[[i]]))
    
    df <- merge(pred, hrdetect.annotations, by='sample', all.x=T)
    df$cohort <- i
    df <- df[order(df['p_chord']),]
    
    return(df)
  }); names(merged_data) <- names(l_pred)
  
  ## Remove BRCA-EU samples not used for performance evaluation
  merged_data$BRCA_EU <- merged_data$BRCA_EU[merged_data$BRCA_EU$used_for_evaluation,]
  
  out <- do.call(rbind,merged_data)
  rownames(out) <- NULL
  
  ## Add features
  sigs_rel2 <- do.call(rbind,l_sigs_rel2)
  
  cbind(
    out,
    sigs_rel2[ match(out$sample,sigs_rel2$sample), ]
  )
}

pred_icgc_pre <- chordPredictIcgc(chord, l_contexts_icgc, hrdetect_annotations)

pred_icgc <- pred_icgc_pre[pred_icgc_pre$cohort!='BRCA_EU.HMF_pipeline',]
pred_icgc.hmf_pipeline <- pred_icgc_pre[pred_icgc_pre$cohort=='BRCA_EU.HMF_pipeline',]

#========= CHORD (sigs model) =========#
#--------- Model ---------#
sig_model <- readRDS(paste0(base_dir,'/HMF_DR010_DR047/scripts/chord_training/main/output/61_snvSigs_svSigs_allSigsRel_noSuppEvidence_noBoruta/seed002/final/rf_out.rds'))$RF

#--------- Outer CV ---------#
outer_cv.sig_model <- gatherOuterCv(paste0(base_dir,'/HMF_DR010_DR047/scripts/chord_training/main/output/61_CV_snvSigs_svSigs_allSigsRel_noSuppEvidence_noBoruta/'))

cv_pred_hmf.sig_model <- gatherOuterCvPred(outer_cv.sig_model)

#--------- HMF predictions ---------#
pred_hmf.sig_model <- chordPredictHmf(sig_model, contexts_hmf, annotation2, mode='sigs')

#--------- ICGC predictions ---------#
pred_icgc.sig_model <- chordPredictIcgc(sig_model, l_contexts_icgc, hrdetect_annotations, mode='sigs')


####################################################################################################
# Export supp tables                                                                               #
####################################################################################################
library(openxlsx)

#========= Predictions and CV output =========#
l_pred_xslx <- list(
  CHORD.HMF = pred_hmf,
  CHORD.HMF.CV = cv_pred_hmf,
  CHORD.BRCA_EU = subset(pred_icgc, cohort=='BRCA_EU' & used_for_evaluation),
  CHORD.BRCA_EU.HMF_pipeline = subset(pred_icgc.hmf_pipeline, cohort=='BRCA_EU.HMF_pipeline' & used_for_evaluation),
  
  sig_model.HMF = pred_hmf.sig_model,
  sig_model.HMF.CV = cv_pred_hmf.sig_model,
  sig_model.BRCA_EU = subset(pred_icgc.sig_model, cohort=='BRCA_EU' & used_for_evaluation)
)


counter <- 0
l_pred_xslx <- lapply(l_pred_xslx, function(i){
  
  df <- i[,c('sample','p_chord','BRCA1','BRCA2','response')]
  df <- df[order(df$p_chord, decreasing=T),]
  
  non_hmf_data_tables <- c('CHORD.BRCA_EU', 'CHORD.BRCA_EU.HMF_pipeline', 'sig_model.BRCA_EU')
  counter <<- counter + 1
  if( !(names(l_pred_xslx)[[counter]] %in% non_hmf_data_tables) ){
    df$sample <- analysis_sample_selection$hmf_id[ match(df$sample, analysis_sample_selection$sample_id) ]
  }

  colnames(df) <- c('sample','prob_hrd','prob_brca1_type_hrd','prob_brca2_type_hrd','brca_deficiency')
  return(df)
})

write.xlsx(l_pred_xslx, paste0(base_dir,'/HMF_DR010_DR047/analysis/chord_perf/supp_tables/pred_all_models_datasets.xlsx'))


#========= Features =========#
contexts2Sigs <- function(contexts, mode='abs'){
  contexts$sv <- contexts$sv[,!grepl('_0e00_',colnames(contexts$sv))]
  transformContexts(
    contexts, simplify.types='indel', lsqnonneg.types=c('snv','sv'), 
    rel.types=if(mode=='abs'){ NULL } else { c('snv','indel','sv') }
  )
}

l_contexts_xslx <- list(
  ## HMF
  HMF.raw = transformContexts(contexts_hmf, simplify.types='indel'),
  HMF.contexts_rel = transformContexts(contexts_hmf, simplify.types=c('snv','indel'), rel.types=c('snv','indel','sv')),
  HMF.sigs_abs = contexts2Sigs(contexts_hmf),
  HMF.sigs_rel = contexts2Sigs(contexts_hmf, mode='rel'),
  
  ## BRCA EU
  BRCA_EU.raw = transformContexts(l_contexts_icgc$BRCA_EU, simplify.types='indel'),
  BRCA_EU.contexts_rel = transformContexts(l_contexts_icgc$BRCA_EU, simplify.types=c('snv','indel'), rel.types=c('snv','indel','sv')),
  BRCA_EU.sigs_abs = contexts2Sigs(l_contexts_icgc$BRCA_EU),
  BRCA_EU.sigs_rel = contexts2Sigs(l_contexts_icgc$BRCA_EU, mode='rel'),
  
  ## BRCA EU HMF pipeline
  BRCA_EU_HMF_pl.raw = transformContexts(l_contexts_icgc$BRCA_EU.HMF_pipeline, simplify.types='indel'),
  BRCA_EU_HMF_pl.contexts_rel = transformContexts(l_contexts_icgc$BRCA_EU.HMF_pipeline, simplify.types=c('snv','indel'), rel.types=c('snv','indel','sv'))
)

counter <- 0
l_contexts_xslx_2 <- lapply(l_contexts_xslx, function(i){

  df <- as.data.frame(i)
  df <- cbind(sample=rownames(df),df)
  rownames(df) <- NULL
  
  counter <<- counter + 1
  table_name <- names(l_contexts_xslx)[[counter]]
  if(grepl('^HMF',table_name)){
    df$sample <- analysis_sample_selection$hmf_id[ match(df$sample, analysis_sample_selection$sample_id) ]
  }
  
  if(grepl('sigs',table_name)){
    colnames(df) <- gsub('^e[.]','cosmic_snv.',colnames(df))
    colnames(df) <- gsub('^SV','RS',colnames(df))
  }
  
  return(df)
})

write.xlsx(l_contexts_xslx_2, paste0(base_dir,'/HMF_DR010_DR047/analysis/chord_perf/supp_tables/features_all_datasets.xlsx'))



####################################################################################################
# Sorted prediction probs                                                                          #
####################################################################################################

plotPredProbs <- function(
  probs, response, sample.names=NA,
  error=NA, error.as.linerange=F,
  
  class.colors=c('BRCA1-type HRD'='#f58225', 'BRCA2-type HRD'='#69439d'),
  response.colors=c('BRCA1'='#f58225', 'BRCA2'='#69439d','none'=NA),
  
  default.class.name='HRD',
  
  cutoff.hline.y=0.5, show.confusion=0, confusion.xpos=NULL,
  rel.heights=c(1, 0.3), ylims=c(0,1), plot.title=NULL, y.lab=NULL, top=NA, bottom=NA, 
  
  hide.info=c(), hide.response=F, do.ordering=T
){
  
  # probs=pred_icgc[c('BRCA1','BRCA2')]
  # colnames(probs) <- c('BRCA1-type HRD','BRCA2-type HRD')
  # response <- as.factor(pred_icgc$response)
  # sample.names=pred_icgc$sample
  
  # error=smooth.spline(cv100_pred_stats$hrd.sd)$y
  # class.colors=c('BRCA1'='#f9b47c', 'BRCA2.mean'='#a58ec4')

  #--------- Pre-process ---------#
  if(!(is.data.frame(probs) || is.matrix(probs))){
    stop('`probs` must be a data.frame or matrix')
  }
  
  if(is.matrix(probs)){ probs <- as.data.frame(probs) }
  df <- cbind(sum=rowSums(probs), probs, response, error, sample.names)
  
  if(do.ordering){
    df <- df[do.call(order, df),]
  }
  
  #df$sum <- NULL
  df$index <- 1:nrow(df)
  
  df_full <- df
  if(!is.na(top)){ 
    df <- tail(df, top) 
  } else if(!is.na(bottom)){ 
    df <- head(df, bottom) 
  }
  
  df_melt <- melt(subset(df, select=-sum),c('index','response','error','sample.names'))
  colnames(df_melt)[5:6] <- c('class','probs')
  
  #--------- Prediction probs ---------#
  plot_pred <- ggplot() +
    geom_bar(
      data=df_melt, mapping=aes(index, probs, fill=class),
      stat='identity', position='stack', width=1
    ) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(limits=ylims) +
    labs(x=paste0(nrow(df_full),' samples'), y='Probability', fill='Prediction class') +
    theme_bw() +
    theme(
      plot.title=element_text(hjust=0.5, size=11),
      axis.title.y=element_text(size=12),
      axis.title.x=element_blank(),
      #axis.ticks.x=element_blank(),
      axis.text.x=element_blank(),
      panel.grid.minor.x=element_blank(),
      panel.grid.major.x=element_blank(),
      panel.grid.minor.y=element_blank(),
      legend.justification=c(0,0.5),
      #legend.key=element_rect(size=1.5, color='white'),
      legend.key.size=unit(11,'pt'),
      legend.text=element_text(size=9),
      legend.title=element_text(size=10, face='bold'),
      plot.margin = unit(c(2,0,6,2),'pt')
    )
  
  if(!is.null(class.colors)){
    plot_pred <- plot_pred + scale_fill_manual(values=class.colors)
  }
  
  if(!anyNA(error)){
    pd_error <- data.frame(
      index=df$index,
      y=rowSums(df[colnames(df) %in% colnames(probs)])
    )
    
    pd_error$upper <- pd_error$y + df$error
    pd_error$lower <- pd_error$y - df$error
    
    if(error.as.linerange){
      plot_pred <- plot_pred + geom_linerange(data=pd_error, aes(index, ymin=lower, ymax=upper), size=0.5, color='black')
    } else {
      pd_error <- melt(pd_error,c('index','y'))
      plot_pred <- plot_pred + geom_line(data=pd_error, aes(index, value, group=variable), size=0.5, color='black')}
  }
  
  if(!is.null(cutoff.hline.y)){
    plot_pred <- plot_pred + geom_hline(yintercept=cutoff.hline.y, linetype='dashed', size=0.3)
  }
  
  if(!is.null(plot.title)){ plot_pred <- plot_pred + ggtitle(plot.title) }
  if(!is.null(y.lab)){ plot_pred <- plot_pred + ylab(y.lab) }
  
  if(show.confusion>=1){
    
    m_confusion <- confusionMatrix(
      df_full$sum, 
      toBinaryResponse(df_full$response,c('BRCA1','BRCA2'),1,'none',0),
      cutoff=cutoff.hline.y
    )
    m_confusion <- rbind(
      m_confusion[c('tp','fp')],
      m_confusion[c('fn','tn')]
    )
    
    rownames(m_confusion) <- NULL
    colnames(m_confusion) <- NULL
    
    if(show.confusion==1){
      cols <- matrix(c('red','grey50'), nrow(m_confusion), ncol(m_confusion), byrow = T)
    } else if(show.confusion==2){
      m_confusion <- as.matrix(m_confusion[,1] + m_confusion[,2])
      cols <- as.matrix(c('#C00000','#548235'))
    }
    
    tt <- ttheme_minimal(core=list(fg_params = list(col = cols)), base_size=10)
    table_confusion <- tableGrob(m_confusion, theme = tt)
    
    if(is.null(confusion.xpos)){
      table_xpos <- nrow(df_full) * 0.1
    } else {
      table_xpos <- confusion.xpos
    }
    
    #plot + annotation_custom(table_confusion, xmin=table_xpos, xmax=table_xpos, ymin=cutoff, ymax=cutoff)
    plot_pred <- plot_pred + 
      annotation_custom(
        table_confusion, xmin=table_xpos, xmax=table_xpos, 
        ymin=cutoff.hline.y, ymax=cutoff.hline.y
      )
  }
  
  #--------- Response ---------#
  df_response <- df_melt[,c('index','response','sample.names')]
  df_response <- df_response[!duplicated(df_response),]
  
  df_response$sample.names <- factor(df_response$sample.names, unique(df_response$sample.names))
  neg_class <- names(response.colors)[is.na(response.colors)]
  levels(df_response$response)[levels(df_response$response) %in% neg_class] <- NA
  # response.colors2 <- response.colors
  # names(response.colors2)[response.colors2 %in% neg_class] <- ''
  
  if(!is.na(sample.names)){
    plot_response <- ggplot(df_response, aes(sample.names, paste0('n = ',nrow(probs)), fill=response))
  } else {
    plot_response <- ggplot(df_response, aes(index, paste0('n = ',nrow(probs)), fill=response)) +
      scale_x_continuous(expand=c(0,0))
  }

  plot_response <- plot_response +
    geom_tile() +
    
    scale_y_discrete(expand=c(0,0)) +
    labs(x=paste0(nrow(df_full),' samples'), fill='Gene deficiency') +
    
    theme_bw() +
    theme(
      axis.line=element_blank(),
      #axis.text.y=element_text(face='bold',size=10),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      axis.title.y=element_blank(),
      axis.title.x=element_text(size=12),
      panel.border=element_rect(fill=NA, color='black', size=0.5),
      panel.grid=element_blank(),
      legend.justification=c(0,0.5),
      legend.key=element_rect(size=1.2, color='white'),
      legend.key.size=unit(10,'pt'),
      legend.text=element_text(size=9),
      legend.title=element_text(size=10, face='bold'),
      plot.margin = unit(c(2,0,2,2),'pt')
    )
  
  if(!is.na(sample.names)){
    plot_response <- plot_response + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
    plot_pred <- plot_pred + theme(axis.ticks.x=element_blank())
  }
  
  if(!is.null(response.colors)){
    plot_response <- plot_response + 
      scale_fill_manual(values=response.colors, breaks=names(response.colors)[!is.na(response.colors)])
  }
  
  #--------- Hide info ---------#
  if('left' %in% hide.info){
    theme_hidden_left <- theme(
      axis.title.y=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()
    )
    
    plot_pred <- plot_pred + theme_hidden_left
    plot_response <- plot_response + theme_hidden_left
  }
  
  if('right' %in% hide.info){
    theme_hidden_right <- theme(legend.position='none')
    plot_pred <- plot_pred + theme_hidden_right
    plot_response <- plot_response + theme_hidden_right
  }
  
  if('bottom' %in% hide.info){
    theme_hidden_bottom <- theme(axis.title.x=element_blank())
    plot_pred <- plot_pred + theme_hidden_bottom
    plot_response <- plot_response + theme_hidden_bottom
  }
  
  #--------- Combine ---------#
  if(hide.response){
    plot_pred + theme(axis.title.x=element_text(), axis.text.x=element_text())
  } else {
    plot_grid(plot_pred, plot_response, align='v', axis='tblr', ncol=1, rel_heights=rel.heights)
  }
}

plotPredProbsHmfSplit <- function(
  pred, top=300, title='CHORD on HMF dataset', 
  show.confusion=0, confusion.xpos=600, cutoff.hline.y=NULL, hide.response=F, axis.title.x=NULL
){
  #pred=cv_pred_hmf
  
  probs <- pred[c('BRCA1','BRCA2')]
  colnames(probs) <- c('BRCA1-type HRD','BRCA2-type HRD')
  
  n_samples <- nrow(pred)
  bottom <- n_samples - top
  
  plots <- list()
  
  plots$left <- plotPredProbs(
    probs, pred$response, plot.title=paste0('Bottom ', bottom), bottom=bottom, 
    hide.info=c('right','bottom'), rel.heights=c(1, 0.25), 
    show.confusion=show.confusion, confusion.xpos=confusion.xpos, cutoff.hline.y=cutoff.hline.y,
    hide.response=hide.response
  )
  
  plots$right <- plotPredProbs(
    probs, pred$response, plot.title=paste0('Top ', top), top=top, 
    hide.info=c('left','bottom'), rel.heights=c(1, 0.25), cutoff.hline.y=cutoff.hline.y,
    hide.response=hide.response
  )
  
  if(hide.response){
    plots <- lapply(plots, function(i){ i + theme(axis.title.x=element_blank()) })
  }
  
  if(is.null(axis.title.x)){
    axis.title.x <- sprintf('%s samples            ', n_samples)
  }
  
  arrangeGrob(plots$left, plots$right, nrow=1, widths=c(0.55, 0.45), top=title, bottom=axis.title.x)
}

#========= Main =========#
##
plotPredProbsWrapper <- function(pred_hmf, cv_pred_hmf, pred_icgc, out_path, model_name='CHORD'){
  
  plots_pred <- list()
  
  #plots_pred$hmf <- plotPredProbsHmfSplit(pred_hmf, title='CHORD on HMF dataset', cutoff.hline.y=0.5, show.confusion=2)
  plots_pred$hmf <- (function(){
    probs <- pred_hmf[pred_hmf$sample %in% selected_samples,]
    plotPredProbsHmfSplit(
      probs, title='CHORD on HMF dataset', cutoff.hline.y=0.5, show.confusion=2, hide.response=T,
      axis.title.x=paste0(nrow(probs),' HMF patients      ')
    )
  })()
  
  plots_pred$hmf_cv <- plotPredProbsHmfSplit(cv_pred_hmf, title='CHORD CV on HMF dataset', confusion.xpos=200)
  
  plots_pred$icgc_chord <- (function(){
    probs <- pred_icgc[c('BRCA1','BRCA2')]
    colnames(probs) <- c('BRCA1-type HRD','BRCA2-type HRD')
    plotPredProbs(probs, pred_icgc$response, plot.title=paste0(model_name,' on ICGC dataset'), cutoff.hline.y=NULL)
  })()
  
  plots_pred$brca_eu_chord <- (function(){
    pred_ss <- pred_icgc[pred_icgc$cohort=='BRCA_EU',]
    probs <- pred_ss[c('BRCA1','BRCA2')]
    response <- pred_ss$response
    colnames(probs) <- c('BRCA1-type HRD','BRCA2-type HRD')
    plotPredProbs(probs, response, plot.title=paste0(model_name,' on BRCA-EU dataset'), cutoff.hline.y=NULL)
  })()
  
  # plots_pred$icgc_chord.hmf_pipeline <- (function(){
  #   probs <- pred_icgc.hmf_pipeline[c('BRCA1','BRCA2')]
  #   colnames(probs) <- c('BRCA1-type HRD','BRCA2-type HRD')
  #   plotPredProbs(probs, pred_icgc.hmf_pipeline$response, plot.title=paste0(model_name,' on ICGC dataset; HMF pipeline'), cutoff.hline.y=NULL)
  # })()
  
  plots_pred$icgc_hrdetect <- (function(){
    probs <- pred_icgc['p_hrdetect']
    colnames(probs) <- 'HRD'
    plotPredProbs(
      probs, pred_icgc$response, plot.title='HRDetect on ICGC dataset',
      class.colors=c(HRD='#bac4ba'), cutoff.hline.y=0.7
    )
  })()
  
  ##
  pdf(out_path, 8, 4)
  counter <- 0
  for(i in plots_pred){ 
    counter <- counter + 1
    if(counter==2){ grid.newpage() } ## Hack to prevent hmf and hmf_cv plots from overlapping on one page
    grid.draw(i)
  }
  dev.off()
  
}

## CHORD
plotPredProbsWrapper(
  pred_hmf, cv_pred_hmf, pred_icgc, 
  out_path=paste0(base_dir,'/HMF_DR010_DR047/analysis/chord_perf/chord_main/sorted_probs2.pdf')
)

## Sigs model
plotPredProbsWrapper(
  pred_hmf.sig_model, cv_pred_hmf.sig_model, pred_icgc.sig_model, 
  paste0(base_dir,'/HMF_DR010_DR047/analysis/chord_perf/sigs_model/sorted_probs.pdf'),
  model_name='Signature model'
)

#========= 100x repeated CV =========#
#--------- mean probs ---------#
plots_cv100 <- (function(){
  l <- list()
  
  ## Probs
  probs <- cv100_pred_stats[c('BRCA1.mean','BRCA2.mean')]
  colnames(probs) <- c('BRCA1-type HRD','BRCA2-type HRD')
  
  error <- smooth.spline(cv100_pred_stats$hrd.sd)$y
  class.colors=c('BRCA1-type HRD'='#f9b47c', 'BRCA2-type HRD'='#a58ec4')
  
  l[[1]] <- plotPredProbs(
    probs, cv100_pred_stats$response, error,
    class.colors=class.colors, ylims=c(NA,NA),
    plot.title='Mean prediction prob. after 100x 10-fold CV on training set',
    y.lab='Mean probability\n(error: SD HRD probability)'
  )
  
  l[[2]] <- plotPredProbs(
    probs, cv100_pred_stats$response, error,
    class.colors=class.colors, ylims=c(NA,NA),
    plot.title='Mean prediction prob. after 100x 10-fold CV on training set (top 200)',
    y.lab='Mean probability\n(error: SD HRD probability)',
    top=200, error.as.linerange=T
  )
  
  return(l)
})()

#--------- n times HRD ---------#
plotTimesHrd <- function(
  cv100_pred_stats,
  response.colors=c('BRCA1'='#f58225', 'BRCA2'='#69439d','none'='grey'),
  
  hline.upper.yintercept=60, hline.upper.label="Below: 'BRCA1' and 'BRCA2' samples filtered",
  hline.lower.yintercept=40, hline.lower.label="Above: 'none' samples filtered",
  
  rel.heights=c(1, 0.2),
  show.top=200
){
  df <- cv100_pred_stats
  df$response <- factor(df$response, c('none','BRCA1','BRCA2'))
  df <- df[do.call(order, df[c('n_is_hrd','response')]),]
  df$index <- 1:nrow(df)
  
  if(!is.na(show.top)){ df <- tail(df, show.top) }
  
  ggplot(df, aes(index, n_is_hrd, color=response)) + 
    geom_point() +
    
    #scale_x_continuous(minor_breaks=df$index) +
    scale_y_continuous(breaks=seq(0, max(df$n_is_hrd), 10)) +
    scale_color_manual(values=response.colors, breaks=names(response.colors)) +
    
    geom_hline(yintercept=hline.upper.yintercept, linetype='dotted') +
    annotate(
      'text', x=min(df$index), y=hline.upper.yintercept, hjust=0, vjust=2, size=2.8,
      label=hline.upper.label
    ) +
    
    geom_hline(yintercept=hline.lower.yintercept, linetype='dotted') +
    annotate(
      'text', x=min(df$index), y=hline.lower.yintercept, hjust=0, vjust=-1, size=2.8,
      label=hline.lower.label
    ) +
    
    labs(
      #y='No. times predicted\nHRD (prob. \u2265 0.5)',
      y=expression("No. times P(HRD) ">=" 0.5"),
      x='Sample', color='Gene deficiency', 
      title='No. times predicted HRD after 100x 10-fold CV on training set (top 200)'
    ) +
    
    theme_bw() +
    theme(
      plot.title=element_text(hjust=0.5, size=11),
      panel.grid.minor.y=element_blank(),
      panel.grid.minor.x=element_blank()
    )
}

plots_cv100[[3]] <- plotTimesHrd(cv100_pred_stats)


pdf(paste0(base_dir,'/HMF_DR010_DR047/analysis/chord_perf/chord_main/sorted_probs_cv100.pdf'), 10, 4)
for(i in plots_cv100){ grid.draw(i) }
dev.off()


####################################################################################################
# Compare CHORD run on samples run through the HMF pipeline vs. downloaded from ICGC               #
####################################################################################################
pred_icgc_origin_vs_hmf_pipeline <- (function(){
  l <- list(
    hmf_pipeline=pred_icgc.hmf_pipeline,
    origin=subset(pred_icgc, cohort=='BRCA_EU')
  )
  
  common_samples <- intersect(l$origin$sample, l$hmf_pipeline$sample)
  
  l <- lapply(l, function(i){
    i[i$sample %in% common_samples,]
  })
  
  l2 <- lapply(l, function(i){
    i[c('sample','BRCA1','BRCA2','p_chord','response')]
  })
  l2$hmf_pipeline <- l2$hmf_pipeline[ order(l2$hmf_pipeline$p_chord), ]
  l2$origin <- l2$origin[ match(l2$hmf_pipeline$sample,l2$origin$sample), ]
  
  ###
  df_melt <- melt(l2, c('sample','p_chord','response'))
  colnames(df_melt)[4:6] <- c('chord_pred_class','prob','cohort')
  
  df_melt <- forceDfOrder(df_melt)
  
  brca_colors <- c(BRCA1='#f58225',BRCA2='#69439d',none=NA)
  
  #========= Sorted probs =========#
  plots <- lapply(unique(df_melt$cohort), function(i){
    
    df_melt_ss <- df_melt[df_melt$cohort==i,]
    
    ggplot(df_melt_ss, aes(sample, prob, fill=chord_pred_class)) +
      #facet_wrap(~cohort, ncol=1) +
      geom_bar(stat='identity', position='stack') +
      scale_fill_manual(name='Prediction class', values=brca_colors, labels=c('BRCA1-type HRD','BRCA2-type HRD')) +
      ylab('Probability') +
      lims(y=c(0,1)) +
      theme_bw() +
      theme(
        panel.grid.minor.y=element_blank(),
        #axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()
      )
  })
  names(plots) <- unique(df_melt$cohort)
  
  # plots$hmf_pipeline <- plots$hmf_pipeline + annotate('text', y=0.75, x=1, vjust=0, hjust=0, label='Variants called from BAM\nfiles with HMF pipeline')
  # plots$origin <- plots$origin + annotate('text', y=0.9, x=1, vjust=0, hjust=0, label='Variants downloaded from ICGC')
  
  plots$hmf_pipeline <- plots$hmf_pipeline + 
    geom_label(aes(x=1, y=1, hjust=0, vjust=1, label='Variants called from BAM\nfiles with HMF pipeline'), fill='white')
  
  plots$origin <- plots$origin +
    geom_label(aes(x=1, y=1, hjust=0, vjust=1, label='Variants downloaded from ICGC'), fill='white')
  
  #========= HRD score difference =========#
  df_diff <- forceDfOrder(data.frame(
    sample=l2$hmf_pipeline$sample,
    diff=l2$hmf_pipeline$p_chord-l2$origin$p_chord
  ))
  
  plots$diff <- ggplot(df_diff, aes(sample,diff)) + 
    geom_bar(stat='identity') +
    geom_hline(yintercept=0, size=0.25, color='red') +
    ylab('Difference') +
    theme_bw() +
    theme(
      panel.grid.minor.y=element_blank(),
      #panel.grid.major.y=element_blank(),
      #axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank()
    )
  
  #========= Response =========#
  df_response <- unique(df_melt[c('sample','response')])
  plots$response <- ggplot(df_response, aes(x=sample,y=1, fill=response)) +
    geom_tile() +
    scale_fill_manual(name='Gene deficiency', values=brca_colors, breaks=c('BRCA1','BRCA2')) +
    scale_x_discrete(expand=c(0,0), name=paste0('\n', nrow(df_response),' ICGC (BRCA-EU) samples')) +
    scale_y_discrete(expand=c(0,0)) +
    theme_bw() +
    theme(
      panel.border=element_rect(fill=NA),
      panel.grid.major=element_blank(),
      axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
      axis.title.y=element_blank()
    )
  
  plots_rearranged <- plots[c('hmf_pipeline','diff','origin','response')]
  
  plot_grid(plotlist=plots_rearranged, align='v',axis='lr',ncol=1, rel_heights=c(1,0.5,1,0.8))
})()

pdf(paste0(base_dir,'/HMF_DR010_DR047/analysis/chord_perf/chord_main/icgc_origin_vs_hmf_pipeline.pdf'), 10, 6)
grid.draw(pred_icgc_origin_vs_hmf_pipeline)
dev.off()


####################################################################################################
# Performance                                                                                      #
####################################################################################################

#========= Confusion =========#
confusionMatrixMcCustom <- function(df, hrd.probs.col='p_chord'){
  confusion <- list()
  
  if(hrd.probs.col=='p_chord'){
    confusion <- confusionMatrixMC(
      df[c('BRCA1','BRCA2','none')], 
      df$response, neg.response='none',cutoff='all'
    )
  }
  
  confusion$hrd <- confusionMatrix( 
    df[[hrd.probs.col]], 
    toBinaryResponse(df$response,c('BRCA1','BRCA2'),1,'none',0) 
  )
  
  confusion <- confusion[c('hrd','BRCA1','BRCA2')]
  names(confusion) <- c('HRD','BRCA1-type HRD','BRCA2-type HRD')
  
  confusion <- confusion[sapply(confusion,function(i){ !is.null(i) })]
  
  return(confusion)
}

confusionMatrixMcWrapper <- function(cv_pred_hmf, pred_icgc){
  l_confusion <- list()
  l_confusion$HMF <- confusionMatrixMcCustom(cv_pred_hmf)
  l_confusion$ICGC_CHORD <- confusionMatrixMcCustom(pred_icgc)
  l_confusion$BRCA_EU_CHORD <- confusionMatrixMcCustom(pred_icgc[pred_icgc$cohort=='BRCA_EU',])
  l_confusion$ICGC_HRDetect <- confusionMatrixMcCustom(pred_icgc, 'p_hrdetect')
  return(l_confusion)
}

#========= Compound performance =========#
plotPerfCompoundCustom <- function(
  confusion, compound.metric='pr', plot.title=NULL,
  #line.colors=c('black','#F8766D','#619CFF'),
  line.colors=c('#7e7e7e','#f58225','#69439d'),
  line.size=c(1, 0.5, 0.5),
  legend.position=c(0.5, 0.1)
){
  
  l <- lapply(confusion, function(i){ 
    m <- calcPerfCompound(i, compound.metric, metric.names.as.x.y=T)
    m <- as.data.frame(m)
    m$auc <- calcAUC(m$x,m$y)
    return(m)
  })
  df <- as.data.frame(do.call(rbind,l))
  rownames(df) <- NULL
  
  ##
  df$class <- unlist(lapply(names(l), function(i){ rep(i, nrow(l[[i]])) }))
  df$class <- sprintf('%s (%s)', df$class, format(round(df$auc,3),nsmall=3))
  df$class <- factor(df$class,unique(df$class))
  
  ##
  axis_titles <- switch(
    compound.metric,
    roc=c('False positive rate','True positive rate'),
    pr=c('Recall (TPR)','Precision (PPV)'),
    npv_tnr=c('True negative rate','Negative predictive value')
  )
  
  plot <- ggplot(df, aes(x,y, color=class, size=class)) + 
    geom_path() +
    labs(x=axis_titles[1],y=axis_titles[2], color='Prediction class (AUC)') +
    theme_bw() +
    theme(
      plot.title=element_text(hjust=0.5, size=10),
      legend.title=element_text(size=10),
      legend.text=element_text(size=9),
      legend.key.height=unit(10,'pt'),
      legend.key.width=unit(12,'pt'),
      legend.position=legend.position,
      legend.justification=c(0.5, 0),
      legend.background=element_rect(color='black',size=0.2,fill=alpha('white',0.8)),
      plot.margin=ggplot2::margin(15,15,15,15)
    )
  
  if(!is.null(line.colors)){ plot <- plot + scale_color_manual(values=line.colors) }
  if(!is.null(line.size)){ plot <- plot + scale_size_manual(values=line.size, guide=F) }
  if(!is.null(plot.title)){ plot <- plot + ggtitle(plot.title) }
  
  return(plot)
}

#========= Simple performance =========#
plotPerfCustom <- function(
  confusion, metrics=c('tpr','fpr'), plot.title=NULL,
  line.colors=NULL, line.types=NULL,
  legend.position=c(0.5, 0.35), na.rough.fix=F
){
  # confusion=l_confusion$HMF
  # metrics=c('tpr','fpr')
  
  l <- lapply(names(confusion), function(i){
    #i='BRCA1-like HRD'
    df <- calcPerf(confusion[[i]], metrics=metrics, melt=T)
    if(na.rough.fix){ df[is.na(df)] <- 0 }
    df <- unique(df)
    df$class <- i
    return(df)
  })
  df <- as.data.frame(do.call(rbind,l))
  rownames(df) <- NULL
  
  ##
  metric_names <- c(
    tpr='True positive rate',
    tnr='True negative rate',
    fpr='False positive rate',
    fnr='False negative'
  )
  
  for(i in names(metric_names)){
    df$metric <- gsub(i, metric_names[i], df$metric, fixed=T)
  }
  
  ##
  for(i in 1:ncol(df)){
    if(is.character(df[[i]])){ df[[i]] <- factor(df[[i]], unique(df[[i]])) }
  }
  
  ##
  plot <- ggplot(df, aes(cutoff, value, color=metric, linetype=class)) + 
    geom_path() +
    labs(x='Classif. cutoff',y='Metric value',color='Metric',linetype='Prediction class') +
    theme_bw() +
    theme(
      plot.title=element_text(hjust=0.5, size=10),
      legend.title=element_text(size=9),
      legend.text=element_text(size=8),
      legend.position=legend.position,
      legend.direction='vertical',
      legend.key.height=unit(10,'pt'),
      legend.key.width=unit(12,'pt'),
      legend.background=element_rect(fill=NA),
      legend.spacing.y=unit(1,'pt'),
      legend.box.background=element_rect(color='black',size=0.2,fill=alpha('white',0.8)),
      plot.margin=ggplot2::margin(15,15,15,15)
    )
  
  if(!is.null(line.colors)){ plot <- plot + scale_color_manual(values=line.colors) }
  if(!is.null(line.types)){ plot <- plot + scale_linetype_manual(values=line.types) }
  if(!is.null(plot.title)){ plot <- plot + ggtitle(plot.title) }
  
  return(plot)
}

#========= Main =========#
plotPerfWrapper <- function(l_confusion, groups=c('HMF','ICGC_CHORD'), out_path){
  
  plots_perf_all <- plot_grid(
    
    ##
    plotPerfCompoundCustom(l_confusion[[ groups[1] ]], 'roc'),
    plotPerfCompoundCustom(l_confusion[[ groups[1] ]], 'pr'),
    plotPerfCustom(l_confusion[[ groups[1] ]]),
    
    ##
    #plotPerfCompoundCustom(l_confusion$ICGC_CHORD, 'roc') +
    plotPerfCompoundCustom(l_confusion[[ groups[2] ]], 'roc') +
      theme(plot.title=element_blank()),
    
    #plotPerfCompoundCustom(l_confusion$ICGC_CHORD, 'pr') +
    plotPerfCompoundCustom(l_confusion[[ groups[2] ]], 'pr') +
      theme(plot.title=element_blank()),
    
    #plotPerfCustom(l_confusion$ICGC_CHORD, legend.position='none') +
    plotPerfCustom(l_confusion[[ groups[2] ]], legend.position='none') +
      theme(plot.title=element_blank()),
    
    nrow=2, ncol=3, align='hv', axis='tblr'
  )
  
  #pdf(paste0(base_dir,'/HMF_DR010_DR047/analysis/chord_perf/chord_main/perf_all.pdf'), 8, 7.5)
  pdf(out_path, 12, 7.5)
  grid.draw(plots_perf_all)
  dev.off()
}

## CHORD
l_confusion <- confusionMatrixMcWrapper(cv_pred_hmf, pred_icgc)
#plotPerfWrapper(l_confusion, out_path=paste0(base_dir,'/HMF_DR010_DR047/analysis/chord_perf/chord_main/perf_all.pdf'))
plotPerfWrapper(l_confusion, groups=c('HMF','BRCA_EU_CHORD'), out_path=paste0(base_dir,'/HMF_DR010_DR047/analysis/chord_perf/chord_main/perf_all_brca_eu_test.pdf'))


## Sigs model
l_confusion.sigs_model <- confusionMatrixMcWrapper(cv_pred_hmf.sig_model, pred_icgc.sig_model)
plotPerfWrapper(l_confusion.sigs_model, out_path=paste0(base_dir,'/HMF_DR010_DR047/analysis/chord_perf/sigs_model/perf_all.pdf'))
plotPerfWrapper(l_confusion.sigs_model, groups=c('HMF','BRCA_EU_CHORD'), out_path=paste0(base_dir,'/HMF_DR010_DR047/analysis/chord_perf/sigs_model/perf_all_brca_eu.pdf'))


####################################################################################################
# Performance on methylation samples                                                               #
####################################################################################################
pred_icgc_meth <- pred_icgc[with(pred_icgc,{ is_somatic_meth & response != 'none' }),]
confusion_meth <- confusionMatrixMcCustom(pred_icgc_meth)

plots_meth <- list()

plots_meth$probs <- (function(){
  probs <- pred_icgc_meth[c('BRCA1','BRCA2')]
  colnames(probs) <- c('BRCA1-type HRD','BRCA2-type HRD')
  plotPredProbs(
    probs, pred_icgc_meth$response, sample.names=pred_icgc_meth$sample,
    plot.title='CHORD on ICGC dataset (samples with methylation)',
    response.colors=c('BRCA1'='#f58225', 'BRCA2'='#69439d','none'=NA), show.confusion=2,
    rel.heights=c(1,0.5)
  )
})()

pdf(paste0(base_dir,'/HMF_DR010_DR047/analysis/chord_perf/chord_main/sorted_probs_meth.pdf'), 8, 4)
grid.draw(plots_meth$probs)
dev.off()


plots_meth$compound <- plotPerfCompoundCustom(confusion_meth)
plots_meth$simple <- plotPerfCustom(confusion_meth, na.rough.fix=T, legend.position='right')

pdf(paste0(base_dir,'/HMF_DR010_DR047/analysis/chord_perf/chord_main/perf_all_meth.pdf'), 9, 3.5)
grid.draw(
  plot_grid(plots_meth$compound, plots_meth$simple, ncol=2, align='v', axis='b', rel_widths=c(1,1.35))
)
dev.off()



####################################################################################################
# Feature importance                                                                               #
####################################################################################################

#========= Functions =========#
getFeatImp <- function(chord, class=NULL, f_feature_rename=NULL){
  df <- as.data.frame(randomForest::importance(chord, class=class, type=1))
  colnames(df) <- 'mda'
  
  feature_names <- rownames(df)
  if(!is.null(f_feature_rename)){
    feature_names <- f_feature_rename(feature_names)
  }
  
  df <- cbind(feature=feature_names,df); rownames(df) <- NULL
  df[order(df$mda, decreasing=F),]
}

getFeatImpCv <- function(chord.cv, class=NULL, melt.df=T, f_feature_rename=NULL){
  #chord.cv=cv_outer_sigs
  
  l <- lapply(chord.cv, function(i){
    #i=chord_cv[[1]]
    chord_fold <- i$RF
    chord_fold$classes
    randomForest::importance(chord_fold, class=class, type=1)[,1]
  })
  
  ## Fill in NA in folds with missing features vs all features across folds
  all_features <- unique(unlist(lapply(l, names), use.names=F))
  l <- lapply(l, function(i){
    #i=l[[1]]
    
    missing_features <- all_features[!(all_features %in% names(i))]
    i[missing_features] <- NA
    
    ## Force consistent order
    return(i[all_features])
  })
  
  df <- as.data.frame(do.call(cbind,l))
  colnames(df) <- 1:length(l)
  df <- df[order(rowMeans(df, na.rm=T), decreasing=F),]
  
  feature_names <- rownames(df)
  if(!is.null(f_feature_rename)){
    feature_names <- f_feature_rename(feature_names)
  }
  
  df <- cbind(feature=feature_names, df); rownames(df) <- NULL
  
  if(!melt.df){
    return(df)
  }
  
  df_melt <- suppressMessages({ melt(df, by='feature') })
  colnames(df_melt)[2:3] <- c('fold','mda')
  return(df_melt)
}

plotFeatImp <- function(imp, imp.cv, plot.title=NULL, hide.axis.title.x=F, seed=1){
  # imp=chord_imp[[1]]
  # imp.cv=chord_cv_imp[[1]]
  
  set.seed(seed)
  
  imp.cv$feature <- factor(imp.cv$feature, unique(imp.cv$feature))
  imp$feature <- factor(imp$feature, unique(imp.cv$feature))
  
  plot <- ggplot(imp.cv, aes(feature, mda)) + 
    geom_boxplot(width=0.8, color='grey', size=0.3, outlier.color=NA) +
    geom_jitter(width=0.1, color='black', size=1, shape=1) +
    geom_point(data=imp, aes(feature, mda), shape=108, size=1, stroke=7, color='red') +
    scale_y_continuous(limits=c(0,NA)) +
    labs(y='Feature importance\n(Mean decrease in accuracy)') +
    
    coord_flip() +
    
    theme_bw() +
    theme(
      plot.title=element_text(hjust=0.5, size=11),
      panel.grid.major.y=element_line(linetype='dashed'),
      axis.title.y=element_blank(),
      axis.title.x=element_text(size=9),
      axis.text.y=element_text(size=9)
    )
  
  if(!is.null(plot.title)){ plot <- plot + ggtitle(plot.title) }
  if(hide.axis.title.x){ plot <- plot + theme(axis.title.x=element_blank()) }
  
  return(plot)
}


#========= Main =========#
plotFeatImpWrapper <- function(model, model_cv, out_dir, pdf.dims=c(11,3), ...){
  rf_classes <- c('HRD','BRCA2','BRCA1')
  l_imp <- lapply(rf_classes, function(i){
    #print(i)
    if(i=='HRD'){
      rf_class <- NULL
    } else {
      rf_class <- i
    }
    
    imp <- list()
    imp$pred <- getFeatImp(model, class=rf_class, ...)
    imp$cv <- getFeatImpCv(model_cv, class=rf_class, ...)
    
    return(imp)
  })
  names(l_imp) <- rf_classes
  
  
  ## 3x1
  plots_imp <- list(
    plotFeatImp(l_imp[[1]]$pred, l_imp[[1]]$cv, plot.title='HRD', hide.axis.title.x=T),
    plotFeatImp(l_imp[[2]]$pred, l_imp[[2]]$cv, plot.title='BRCA2-type HRD'),
    plotFeatImp(l_imp[[3]]$pred, l_imp[[3]]$cv, plot.title='BRCA1-type HRD', hide.axis.title.x=T)
  )
  plots_imp_combined <- plot_grid(plotlist=plots_imp, nrow=1, align='hv', axis='tblr')
  
  #pdf(paste0(base_dir,'/HMF_DR010_DR047/analysis/chord_perf/chord_main/feat_imp.pdf'), 11, 3)
  pdf(paste0(out_dir,'/feat_imp.pdf'), pdf.dims[1], pdf.dims[2])
  grid.draw(plots_imp_combined)
  dev.off()
}

## CHORD
renameFeaturesContexts <- function(x){
  which_base_subs <- grep('^[A-Z][.][A-Z]',x)
  x[which_base_subs] <- gsub('[.]','>',x[which_base_subs])
  return(x)
}

plotFeatImpWrapper(
  chord, outer_cv, 
  out_dir=paste0(base_dir,'/HMF_DR010_DR047/analysis/chord_perf/chord_main/'),
  f_feature_rename=renameFeaturesContexts
)

## Sigs model ##
renameFeaturesSigs <- function(x){
  which_base_subs <- grep('^e[.]',x)
  x[which_base_subs] <- gsub('^e','cosmic_snv',x[which_base_subs])
  
  # which_sv <- grep('^SV', x)
  # x[which_sv] <- gsub('SV','RS',x[which_sv])
  
  return(x)
}


plotFeatImpWrapper(
  sig_model, outer_cv.sig_model, 
  out_dir=paste0(base_dir,'/HMF_DR010_DR047/analysis/chord_perf/sigs_model/'),
  f_feature_rename=renameFeaturesSigs,
  pdf.dims=c(11,4)
)
