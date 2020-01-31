library(mutSigExtractor)
library(reshape2)
library(grid)
library(gridExtra)
library(ggplot2)
library(cowplot); theme_set(theme_grey())

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
# Data prep                                                                                        #
####################################################################################################

#========= Load data =========#
diplotypes <- read.delim(paste0(base_dir,'/HMF_DR010_DR047/analysis/pancancer_overview/data/diplotypes_filt.txt.gz'))

## Rm duplicate biopsies and keep samples with >=50 indels
analysis_sample_selection <- read.delim(paste0(base_dir,'/HMF_DR010_DR047/analysis/pancancer_overview/data/analysis_sample_selection.txt'))
selected_samples <- analysis_sample_selection[analysis_sample_selection$is_selected,'sample_id']

# PCa_samples <- metadata[metadata$primary_tumor_location=='Prostate','sample_id']

## Based on: https://doi.org/10.1016/j.ygyno.2015.02.017
#hr_panel_genes <- c('BRCA1','BRCA2','ATM','BARD1','BRIP1','CHEK2','MRE11A','NBN','PALB2','RAD50','RAD51C','RAD51D','XRCC2')
hr_panel_genes <- c('BRCA1','BRCA2','RAD51C','PALB2')

#diplotypes <- diplotypes[diplotypes$hgnc_symbol %in% hr_panel_genes,]
#diplotypes$hgnc_symbol <- factor(diplotypes$hgnc_symbol, levels=hr_panel_genes)

#eff_metadata <- read.delim(paste0(base_dir,'/HMF_DR010_DR047/analysis/pancancer_overview/scripts/eff_metadata.txt'))

pred <- (function(){
  df <- read.delim(paste0(base_dir,'/HMF_DR010_DR047/scripts/chord_training/main/output/60_snvContext_svContext_allSigsRel_noSuppEvidence_customBlacklist3_noBoruta/seed015/whole_hmf_dataset_probs.txt'))
  df <- cbind(sample=rownames(df),df)
  rownames(df) <- NULL
  df <- df[df$sample %in% selected_samples,]
  return(df)
})()

#--------- Constants ---------#
hrd_samples <- pred$sample[pred$hrd>=0.5]
hrp_samples <- pred$sample[pred$hrd<0.5]

#chord_group_counts <- list(hrd=length(hrd_samples), hrp=length(hrp_samples))



#========= Determine panel outcomes =========#
#--------- Main ---------#
summarizeDiplotypes <- function(diplotypes, genes=c('BRCA1','BRCA2','RAD51C','PALB2')){
  df <- diplotypes[
    diplotypes$hgnc_symbol %in% genes & diplotypes$diplotype_origin != 'germ_som'
  ,]
  
  df <- df[
    df$hgnc_symbol %in% genes,
    c('sample','hgnc_symbol','diplotype_origin','a1','a2','a2.max_score','a2.max_score_origin')
  ]
  df$a2.origin <- sapply(strsplit(as.character(df$diplotype_origin),'_'),`[`,2, USE.NAMES=F) ## Allele a2 has all germline/somatic mutations
  
  ##
  has_lp_p <- as.integer(
    (df$a2.origin!='cnv' & df$a2.max_score==5) ## snpeff predicted frameshift
    | ( df$a2.max_score==4 & grepl('clinvar',df$a2.max_score_origin) )
  )
  
  ##
  has_lp_p.germ <- as.integer(has_lp_p & df$a2.origin=='germ')
  has_lp_p.som <- as.integer(has_lp_p & df$a2.origin=='som')
  
  has_loh <- as.integer(df$a1=='loh')
  
  ##
  is_biall <- as.integer(
    (df$a1=='loh' & has_lp_p)
    | df$a1 %in% c('full_gene_loss','trunc')
  )
  
  out <- df[,c('sample','hgnc_symbol')]
  out <- cbind(out, has_lp_p.germ, has_lp_p.som, has_loh, is_biall)
  
  ## Ensure all samples have an entry for each gene
  l <- split(out,out$sample)
  out <- do.call(rbind, lapply(l, function(i){
    i <- i[1:4,]
    i$sample <- i$sample[1]
    i$hgnc_symbol[is.na(i$hgnc_symbol)] <- genes[!(genes %in% i$hgnc_symbol)]
    return(i)
  }))
  
  out[is.na(out)] <- 0
  rownames(out) <- NULL
  
  return(out)
}

diplotypes_summary <- summarizeDiplotypes(diplotypes)

#--------- Add metadata ---------#
diplotypes_summary$is_hrd <- as.integer(diplotypes_summary$sample %in% hrd_samples)
diplotypes_summary$cancer_type <- analysis_sample_selection$primary_tumor_location[ match(diplotypes_summary$sample, analysis_sample_selection$sample_id) ]

#--------- Non-cumulative ---------#
detPanelOutcomes <- function(
  diplotypes_summary, 
  group_totals=c(hrd=length(hrd_samples), hrp=length(hrp_samples))
){
  
  is_hrd <- unique(diplotypes_summary[,c('sample','is_hrd')])$is_hrd
  
  if(is.null(group_totals)){
    group_totals=c(
      hrd=sum(is_hrd==1), 
      hrp=sum(is_hrd==0)
    )
  }
  
  panel_positives <- with(diplotypes_summary,{
    list(
      brca.lp_p_germ =
        hgnc_symbol %in% c('BRCA1','BRCA2')
        & has_lp_p.germ
      ,
      
      brca.lp_p_germ_som = 
        hgnc_symbol %in% c('BRCA1','BRCA2') 
        & (has_lp_p.germ | has_lp_p.som) 
      ,
      
      brca.biall = 
        hgnc_symbol %in% c('BRCA1','BRCA2') 
        & is_biall
      ,
      
      hr.biall = is_biall==1
    )
  })
  
  l <- list()
  
  l$hrd <- (function(){
    df <- data.frame(
      neg=NA,
      pos=sapply(panel_positives, function(i){ sum(i & diplotypes_summary$is_hrd==1) })
    )
    df$neg <- group_totals['hrd'] - df$pos
    return(df)
  })()
  
  l$hrp <- (function(){
    df <- data.frame(
      neg=NA,
      pos=sapply(panel_positives, function(i){ sum(i & diplotypes_summary$is_hrd==0) })
    )
    df$neg <- group_totals['hrp'] - df$pos
    return(df)
  })()
  
  return(l)
  
}

l_panel_outcomes <- list(
  BrOv=detPanelOutcomes(subset(diplotypes_summary, cancer_type %in% c('Breast','Ovary'))),
  pancancer=detPanelOutcomes(diplotypes_summary)
)

#--------- Cumulative ---------#
detPanelOutcomesCum <- function(
  diplotypes_summary,
  sel_cancer_types=list(c('Breast','Ovary')),
  group_totals=c(hrd=length(hrd_samples), hrp=length(hrp_samples))
){
  is_hrd <- unique(diplotypes_summary[,c('sample','is_hrd')])$is_hrd
  
  if(is.null(group_totals)){
    group_totals=c(
      hrd=sum(is_hrd==1), 
      hrp=sum(is_hrd==0)
    )
  }
  
  getPanelPositives <- function(diplotypes_summary, infer_group_totals=T){
    
    with(diplotypes_summary,{
      
      df <- data.frame(
        brca.lp_p_germ =
          hgnc_symbol %in% c('BRCA1','BRCA2')
          & has_lp_p.germ
        ,
        
        brca.lp_p_germ_som = 
          hgnc_symbol %in% c('BRCA1','BRCA2') 
          & (has_lp_p.germ | has_lp_p.som) 
        ,
        
        brca.biall = 
          hgnc_symbol %in% c('BRCA1','BRCA2') 
          & is_biall
        ,
        
        hr.biall = is_biall==1
      )
      
      pos <- data.frame(
        hrd=colSums(df & is_hrd==1),
        hrp=colSums(df & is_hrd==0)
      )
      
      if(infer_group_totals){
        df_is_hrd <- unique(diplotypes_summary[,c('sample','is_hrd')])
        group_totals <- c(
          hrd=sum(df_is_hrd$is_hrd==1), 
          hrp=sum(df_is_hrd$is_hrd==0)
        )
      }
      
      neg <- data.frame(
        hrd=group_totals['hrd'] - pos$hrd,
        hrp=group_totals['hrp'] - pos$hrp,
        row.names=rownames(pos)
      )
      
      list(
        pos=as.matrix(pos),
        neg=as.matrix(neg)
      )
    })
    
  }
  #getPanelPositives(diplotypes_summary)
  #subset(diplotypes_summary, hgnc_symbol %in% c('BRCA1','BRCA2') & is_biall & is_hrd)
  #subset(diplotypes_summary, cancer_type %in% c('Breast','Ovary') & hgnc_symbol %in% c('BRCA1','BRCA2') & has_lp_p.germ & is_hrd==0)
  
  remain_cancer_types <- unique(diplotypes_summary$cancer_type)
  
  ## Selected cancer types
  l <- list()
  for(i in sel_cancer_types){
    df <- subset(diplotypes_summary, cancer_type %in% i)
    remain_cancer_types <- remain_cancer_types[!(remain_cancer_types %in% i)]
    
    i_name <- paste(i, collapse='_')
    
    l[[i_name]] <- getPanelPositives(df)
  }
  
  
  
  ## Remaining cancer types
  if( length(remain_cancer_types)>0 ){
    df <- subset(diplotypes_summary, cancer_type %in% remain_cancer_types)
    l[['remainder']] <- getPanelPositives(df)
  }
  
  return(l)
}

l_panel_outcomes_cum <- detPanelOutcomesCum(diplotypes_summary)

####################################################################################################
# Plotting                                                                                         #
####################################################################################################
THEME_DEFAULT <- 
  theme_grey() +
  #theme_classic() +
  theme(
    panel.border=element_rect(fill=NA, color='grey'),
    #panel.border=element_blank(),
    panel.grid.minor.y=element_blank(),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.x=element_blank(),
    
    #legend.title=element_blank(),
    #legend.spacing.x=unit(4,'pt'),
    legend.key=element_rect(color=NA),
    #legend.key.width=unit(1,'lines'),
    #legend.key.height=unit(2,'lines'),
    legend.justification=c(0,0.5),
    
    axis.text.x=element_text(angle=25, hjust=1, vjust=1),
    axis.title.x=element_blank()
  )

#========= % patients missed by genetic screening =========#
plot_panel_outcomes <- (function(){
  
  pd <- (function(){
    pancancer <- l_panel_outcomes$pancancer$hrd[,'pos']
    BrOv <- l_panel_outcomes$BrOv$hrd[,'pos']
    
    pancancer_gain <- pancancer - BrOv
    remainder <- length(hrd_samples) - (BrOv + pancancer_gain)
    
    df <- data.frame(
      remainder=remainder,
      pancancer_gain=pancancer_gain,
      BrOv=BrOv
    )
    rownames(df) <- rownames(l_panel_outcomes$BrOv$hrd)
    
    return(df)
  })()
  
  pd_melt <- melt(as.matrix(pd))
  colnames(pd_melt) <- c('panel_type','bar_fill','abs')
  pd_melt$rel <- melt(pd/length(hrd_samples))$value
  
  pd_melt$labels <- with(pd_melt,{
    paste0(
      signif(100*rel,3),'%',
      '\n(',abs,')'
    )
  })
  
  ggplot(pd_melt, aes(x=panel_type,y=rel, fill=bar_fill, label=labels)) +
    #facet_wrap(cancer_type~., nrow=2) +
    geom_bar(stat='identity', position=position_stack(), color='black', size=0.25) +
    geom_text(size=2.7, position=position_stack(vjust=0.5)) +
    
    scale_y_continuous(
      name='% of HRD tumors',
      labels=function(x){ paste0(x*100,'%') }
    ) +
    
    scale_fill_manual(
      values=c('#B1C4E4','#ede4d8','#F8D598'),
      name=c('Pathogenic event found\nby genetic screening?'),
      labels=c(
        'No',
        'Yes, in other tumor types',
        'Yes, in breast/ovarian tumors'
      )
    ) +
    
    guides(fill=guide_legend(reverse=T)) +
    
    THEME_DEFAULT
})()

plot_panel_outcomes2 <- (function(){
  
  l <- list()
  
  l$abs <- lapply(l_panel_outcomes,`[[`,'hrd')
  l$rel <-   lapply(l$abs, function(i){ i / length(hrd_samples) })
  
  pd <- lapply(l, function(i){
    i <- lapply(i, function(j){
      j <- cbind(panel_type=rownames(j), j)
      rownames(j) <- NULL
      return(j)
    })
    
    i <- do.call(rbind, i)
    i <- cbind(
      cancer_type=sapply(strsplit(rownames(i),'[.]'),`[[`,1),
      i
    )
    rownames(i) <- NULL
    return(i)
  })
  
  pd_melt_pre <- lapply(pd, function(i){
    #id_cols <- c('cancer_type','panel_type')
    i_melt <- melt(i, c('cancer_type','panel_type'))
    i_melt <- forceDfOrder(i_melt)
    i_melt[order(i_melt$cancer_type, i_melt$panel_type),]
  })
  
  pd_melt <- pd_melt_pre$abs
  colnames(pd_melt)[3:4]  <- c('panel_outcome','abs')
  pd_melt$rel <- pd_melt_pre$rel$value
  
  pd_melt$labels <- with(pd_melt,{
    paste0(
      signif(100*rel,3),'%',
      '\n(',abs,')'
    )
  })
  
  pd_melt$fill_factor <- with(pd_melt,{
    paste0(cancer_type, ':', panel_outcome)
  })
  
  fill_pars <- data.frame(
    id=c('BrOv:neg','BrOv:pos','pancancer:neg','pancancer:pos'),
    color=c('#ede4d8','#F8D598', '#cccae0', '#9e9ac4'),
    label=c('No','Yes','   ','   '),
    stringsAsFactors=F
  )

  ggplot(pd_melt, aes(x=cancer_type, y=rel, fill=fill_factor, label=labels)) +
    facet_wrap(.~panel_type, nrow=1, strip.position='bottom') +
    geom_bar(stat='identity', color='black', size=0.25, width=0.9) +
    geom_text(position=position_stack(vjust=0.5), size=3) +
    
    scale_y_continuous(
      name='% of HRD tumors',
      labels=function(x){ paste0(x*100,'%') }
    ) +
    scale_x_discrete(expand=c(0.4,0.4)) +
    
    scale_fill_manual(
      name='Pathogenic event found\nby genetic screening?\n\n        BrOv   Pancancer',
      values=structure(fill_pars$color, names=fill_pars$id),
      labels=fill_pars$label
    ) +
    
    guides(
      fill=guide_legend(
        ncol=2, label.position='left', title.theme=element_text(size=9),
        keywidth=1.5, keyheight=1.5
      )
    ) +
    
    THEME_DEFAULT +
    theme(
      panel.spacing.x=unit(-0.2,'pt'),
      #panel.border=element_blank(),
      strip.placement='outside',
      strip.background=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank()
    )
})()

plot_panel_outcomes3 <- (function(){
  
  #--------- Prep data ---------#
  df_melt <- melt(l_panel_outcomes_cum)
  colnames(df_melt) <- c('panel_type','hr_status','abs','panel_outcome','cancer_type')
  
  #--------- Make labels ---------#
  df_melt <- within(df_melt,{
    
    rel <- NA
    
    rel[hr_status=='hrd'] <- abs[hr_status=='hrd'] / length(hrd_samples)
    rel[hr_status=='hrp'] <- abs[hr_status=='hrp'] / length(hrp_samples)
  })
  
  df_melt <- within(df_melt,{
    perc <- rel*100
    perc <- ifelse(perc>1, round(perc,1), round(perc,2))
  })

  df_melt$labels <- with(df_melt,{
    paste0(
      abs,
      ifelse(
        perc>1,
        paste0('\n(',perc,'%)'),
        paste0(' (',perc,'%)')
      )
    )
  })
  
  #--------- Set ordering ---------#
  df_melt <- df_melt[with(df_melt,{ order(hr_status, cancer_type) }), ]
  df_melt <- df_melt[with(df_melt,{ !(hr_status=='hrp' & panel_outcome=='neg') }), ]
  df_melt[df_melt$hr_status=='hrp','abs'] <- -1 * df_melt[df_melt$hr_status=='hrp','abs']
  
  df_melt$fill_factor <- with(df_melt,{
    paste0(cancer_type,'.',panel_outcome)
  })
  
  df_melt$order_factor <- with(df_melt,{
    paste0(hr_status,'.',fill_factor)
  })
  
  df_melt <- do.call(rbind,split(df_melt,df_melt$order_factor)[
    c(
      'hrd.remainder.pos','hrd.remainder.neg','hrd.Breast_Ovary.neg','hrd.Breast_Ovary.pos',
      #'hrd.remainder.neg','hrd.Breast_Ovary.neg','hrd.remainder.pos','hrd.Breast_Ovary.pos',
      'hrp.Breast_Ovary.pos','hrp.remainder.pos'
    )
  ])
  
  #df_melt$facet_labels <- paste0('CHORD-', toupper(df_melt$hr_status))
  df_melt$facet_labels <- toupper(df_melt$hr_status)
  
  df_melt <- forceDfOrder(df_melt)
  
  #--------- Plot ---------#
  fill_pars <- data.frame(
    name=c('Breast_Ovary.pos','Breast_Ovary.neg','remainder.pos','remainder.neg'),
    color=c('#F8D598', '#f4eee7', '#9e9ac4', '#e0dfec'),
    labels=c('Yes','No','   ','   ')
  )
  
  p <- ggplot(df_melt, aes(panel_type, abs, fill=fill_factor, label=labels)) +
    facet_grid(facet_labels~., scales='free_y', space='free_y', switch='y') +
    #facet_wrap(facet_labels~., scales='free_y', ncol=1) +
    
    geom_bar(stat='identity', color='grey50', size=0.25) +
    geom_text(position=position_stack(vjust=0.5), size=2.5) +
    #geom_text_repel(position=position_stack(vjust=0.5), size=2.5, direction='y', seed=1) + 
  
    scale_y_continuous(
      name='Number of patientss\n(% of patients within group)', expand=c(0,7), 
      labels=function(x){ abs(x) },
      breaks=function(x){ if(max(x)>30){ seq(0,length(hrd_samples),50) } else { c(-40,-20,0) } } ## Hack to get different breaks for HRD and HRP facets
    ) +
    
    scale_fill_manual(
      name='Pathogenic event found\nby genetic testing?',
      values=structure(fill_pars$color, names=fill_pars$name),
      breaks=fill_pars$name,
      labels=fill_pars$labels
    ) +
    
    guides(
      fill=guide_legend(
        ncol=2, label.position='left', title.theme=element_text(size=10),
        keywidth=1.5, keyheight=1.5
      )
    ) +
    
    THEME_DEFAULT +
    theme(
      #strip.placement='outside',
      #strip.background=element_blank(),
      strip.background=element_rect(color='black',size=0.25),
      strip.text=element_text(size=10),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      #axis.title.y=element_blank()
    )
  
  plot_code <- ggplot_gtable(ggplot_build(p))
  
  #--------- Change strip colors ---------#
  strip_name <- which(grepl('strip-l', plot_code$layout$name))
  strip_fills <- c('#F0B7B0','#BCE1CE')
  
  k <- 1
  for (i in strip_name) {
    j <- which(grepl('rect', plot_code$grobs[[i]]$grobs[[1]]$childrenOrder))
    plot_code$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- strip_fills[k]
    k <- k+1
  }
  
  plot_code
  
})()

#========= FDR =========#
plotPanelFdr <- function(bar.width=0.9){
  calcFdr <- function(panel_outcomes){
    #panel_outcomes=l_panel_outcomes$pancancer
    tp <- panel_outcomes$hrd$pos
    fp <- panel_outcomes$hrp$pos   
    
    #panel_outcomes <- l_panel_outcomes_cum$Breast_Ovary
    # tp <- panel_outcomes$pos[,'hrd']
    # fp <- panel_outcomes$pos[,'hrp']
    
    df <- data.frame(
      panel_type=rownames(panel_outcomes$hrd), 
      #panel_type=rownames(panel_outcomes$pos), 
      tp, 
      fp, 
      total=tp+fp
    )
    
    df$fdr <- df$fp/df$total
    
    df$label <- with(df,{
      paste0(
        signif(100*fdr,3),'%',
        '\n(',fp,'/',total,')'
      )
    })
    
    return(df)
  }

  pd <- do.call(rbind,lapply(l_panel_outcomes, calcFdr))
  pd <- cbind(
    cancer_type=sapply(strsplit(rownames(pd),'[.]'),`[[`,1),
    pd
  )
  rownames(pd) <- NULL
  
  pd <- forceDfOrder(pd)
  
  ggplot(pd, aes(x=panel_type, y=fdr, fill=cancer_type, label=label)) +
    #facet_wrap(.~cancer_type, ncol=1) +
    geom_bar(
      stat='identity', position=position_dodge(width=bar.width), width=bar.width, 
      color='black', size=0.25
    ) +
    geom_text(position=position_dodge(width=bar.width), size=2.7) +
    
    scale_y_continuous(
      name='% of tumors with pathogenic\nevent that are CHORD-HRP',
      labels=function(x){ paste0(x*100,'%') },
      expand=c(0,0.015,0,0.025)
    ) +
    
    scale_fill_manual(
      name='Tumor types\nscreened',
      #values=c('#F8D598','#9e9ac4'), 
      values=c('#F8D598','#beafb5'), 
      labels=c('Breast/ovarian','Pan-cancer')
    ) +
    
    THEME_DEFAULT +
    theme(
      axis.ticks.x=element_blank(),
      axis.text.x=element_blank()
    )
  
}

#========= Export =========#
# pdf(paste0(base_dir,'/HMF_DR010_DR047/analysis/pancancer_overview/plots/panel_outcomes.pdf'), 7.5, 4*2)
# grid.draw(
#   plot_grid(
#     plot_panel_outcomes, plotPanelFdr(),
#     ncol=1, align='v'
#   )
# )
# dev.off()


# pdf(paste0(base_dir,'/HMF_DR010_DR047/analysis/pancancer_overview/plots/panel_outcomes2.pdf'), 9, 4*2)
# grid.draw(
#   plot_grid(
#     plot_panel_outcomes2, plotPanelFdr(bar.width=0.75),
#     ncol=1, align='v', axis='tblr'
#   )
# )
# dev.off()

####
# pdf(paste0(base_dir,'/HMF_DR010_DR047/analysis/pancancer_overview/plots/panel_outcomes3.pdf'), 8, 4*2)
# plot_grid(
#   plot_panel_outcomes3, plotPanelFdr(bar.width=0.75),
#   ncol=1, align='v', axis='tblr', rel_heights=c(0.6, 0.4)
# )
# dev.off()

####
pdf(paste0(base_dir,'/HMF_DR010_DR047/analysis/pancancer_overview/plots/panel_outcomes3.pdf'), 8, 7)
plot_grid(
  plot_panel_outcomes3, plotPanelFdr(bar.width=0.75),
  ncol=1, align='v', axis='tblr', rel_heights=c(1, 0.5)
)
dev.off()


# pdf(paste0(base_dir,'/HMF_DR010_DR047/analysis/pancancer_overview/plots/panel_outcomes3_short.pdf'), 8, 4)
# plot(plot_panel_outcomes3)
# dev.off()


























# 
# 
# panel_outcomes <- detPanelOutcomes(diplotypes_summary)
# calcRelPanelOutcomes <- function(
#   panel_outcomes,
#   group_totals=c(hrd=length(hrd_samples), hrp=length(hrp_samples))
# ){
#   
# }
# 
# lapply(l_panel_outcomes, )
# 
# #========= Plot % patients missed by genetic screening =========#
# meltPanelOutcomes <- function(panel_outcomes, mk.labels=T){
#   
#   suppressMessages({
#     l_melt <- list(
#       abs=melt(lapply(panel_outcomes$abs, as.matrix)),
#       rel=melt(lapply(panel_outcomes$rel, as.matrix))
#     )
#   })
#   
#   l_melt <- lapply(l_melt, function(i){
#     #i <- l_melt[[1]]
#     colnames(i) <- c('panel_type','panel_outcome','value','hr_state')
#     i <- i[c(4,1,2,3)]
#     i[order(i$hr_state, i$panel_type),]
#   })
#   
#   df <- l_melt$abs
#   
#   colnames(df)[colnames(df)=='value'] <- 'abs'
#   df$rel <- l_melt$rel$value
#   
#   if(mk.labels){
#     df$labels <- with(df,{
#       paste0(
#         signif(100*rel,3),'%',
#         '\n(',abs,')'
#       )
#     })
#   }
#   
#   return(df)  
# }
# 
# panel_outcomes_merged <- do.call(rbind,lapply(l_panel_outcomes, meltPanelOutcomes))
# 
# panel_outcomes_merged <- cbind(
#   cancer_type=sapply(strsplit(rownames(panel_outcomes_merged),'[.]'),`[[`,1),
#   panel_outcomes_merged
# )
# rownames(panel_outcomes_merged) <- NULL
# 
# panel_outcomes_merged <- as.data.frame(lapply(panel_outcomes_merged, function(i){
#   if(!is.numeric(i)){ i <- factor(i, unique(i)) }
#   return(i)
# }))
# 
# pd <- subset(panel_outcomes_merged, hr_state=='hrd')
# 
# ggplot(pd, aes(x=panel_type,y=rel, fill=panel_outcome, group=cancer_type, label=labels)) +
#   facet_wrap(cancer_type~., nrow=2) +
#   geom_bar(stat='identity', position=position_stack(), color='black', size=0.25) +
#   geom_text(size=2.7, position=position_stack(vjust=0.5)) +
#   scale_y_reverse() +
#   scale_fill_manual(values=c('#ede4d8','#F8D598')) +
#     theme_grey() +
#     theme(
#       panel.border=element_rect(fill=NA, color='black'),
#       panel.grid.minor.y=element_blank(),
#       panel.grid.minor.x=element_blank(),
#       panel.grid.major.x=element_blank(),
# 
#       #legend.title=element_blank(),
#       #legend.spacing.x=unit(4,'pt'),
#       legend.key=element_rect(color=NA),
#       #legend.key.width=unit(1,'lines'),
#       #legend.key.height=unit(2,'lines'),
#       legend.justification=c(0,0.5),
# 
#       axis.text.x=element_text(angle=25, hjust=1, vjust=1),
#       axis.title.x=element_blank()
#     )
# 
# 
# # pd <- subset(panel_outcomes_merged, hr_state=='hrd' & cancer_type=='BrOv')
# # ggplot(pd, aes(x=panel_type,y=rel, fill=panel_outcome)) + 
# #   geom_bar(stat='identity', color='black', size=0.25) +
# #   scale_fill_manual(values=c('#ede4d8','#F8D598')) +
# #   labs()
# 
# # pd <- subset(panel_outcomes_merged, hr_state=='hrd' & cancer_type=='pancancer')
# # ggplot(pd, aes(x=panel_type,y=rel, fill=panel_outcome)) + 
# #   geom_bar(stat='identity')
# 
# 
# 
# 
# 
# panel_outcomes <- detPanelOutcomes(diplotypes_summary)
# 
# 
# # tab <- table(metadata$primary_tumor_location)
# # sel_cancer_types <- c('Breast','Ovary')
# # sum(tab[sel_cancer_types])
# # sum(tab[!(names(tab) %in% sel_cancer_types)])









# panel_outcomes <- with(diplotypes_summary,{
#   # n_samples <- length(unique(sample))
#   # n_hrd <- length(unique(sample[is_hrd==1]))
#   # n_hrp <- length(unique(sample[is_hrd==0]))
#   
#   panel_positives <- list(
#     BrOv.brca_lp_p_germ =
#       cancer_type %in% c('Breast','Ovarian')
#       & hgnc_symbol %in% c('BRCA1','BRCA2')
#       & has_lp_p.germ
#     ,
#     
#     BrOv.brca_lp_p_germ_som = 
#       cancer_type %in% c('Breast','Ovarian')
#       & hgnc_symbol %in% c('BRCA1','BRCA2') 
#       & (has_lp_p.germ | has_lp_p.som) 
#     ,
#     
#     # BrOv.brca_lp_p_germ_som_loh = 
#     #   cancer_type %in% c('Breast','Ovarian')
#     #   & hgnc_symbol %in% c('BRCA1','BRCA2') 
#     #   & (has_lp_p.germ | has_lp_p.som | has_loh) 
#     # ,
#     
#     # BrOv.brca_biall =
#     #   cancer_type %in% c('Breast','Ovarian')
#     #   & hgnc_symbol %in% c('BRCA1','BRCA2')
#     #   & is_biall
#     # ,
#     
#     pancancer.brca_lp_p_germ_som = 
#       hgnc_symbol %in% c('BRCA1','BRCA2') 
#       & (has_lp_p.germ | has_lp_p.som) 
#     ,
#     
#     pancancer.brca_biall = 
#       hgnc_symbol %in% c('BRCA1','BRCA2') 
#       & is_biall
#     ,
#     
#     pancancer.hr_genes_biall = is_biall==1
#   )
#   
#   data.frame(
#     hrd=sapply(panel_positives, function(i){ sum(i & is_hrd==1) }),
#     hrp=sapply(panel_positives, function(i){ sum(i & is_hrd==0) })
#   )
# })
# 
# #========= Plotting =========#
# plot_panel_confusion <- (function(){
#   
#   #panel_outcomes
#   group_totals <- list(
#     hrd=length(pred$sample[pred$hrd>=0.5]),
#     hrp=length(pred$sample[pred$hrd<0.5])
#   )
#   
#   # if(names(panel_outcomes) != names(group_totals)){
#   #   stop("names(panel_outcomes) != names(group_totals)")
#   # }
#   
#   ## Calculate % per HRD/HRP group
#   l <- lapply(names(panel_outcomes), function(i){
#     panel_positive <- panel_outcomes[[i]]
#     panel_negative <- group_totals[[i]] - panel_positive
#     
#     abs <- cbind(panel_negative, panel_positive)
#     rel <- abs/rowSums(abs)
#     
#     list(abs=abs, rel=rel)
#   })
#   
#   names(l) <- names(panel_outcomes)
#   
#   ## Initialize plot data
#   df_melt <- melt(lapply(l,`[[`,'abs'))
#   colnames(df_melt) <- c('panel_type','panel_outcome','abs','is_hrd')
#   
#   df_melt$group <- paste0(df_melt$is_hrd,'.',df_melt$panel_outcome)
#   df_melt$index <- 1:nrow(df_melt)
#   
#   ## Make bar labels
#   df_melt$rel <- melt(lapply(l,`[[`,'rel'))$value
#   
#   df_melt$bar_labels <- with(df_melt,{
#     paste0(
#       abs,
#       '\n(',signif(100*rel,3),'%',')'
#     )
#   })
#   
#   ## Remove panel true negatives
#   df_melt <- subset(df_melt, !(is_hrd=='hrp' & panel_outcome=='panel_negative'))
#   
#   ## Replace % counts with false discovery rate (FDR)
#   fp_sel <- with(df_melt,{ is_hrd=='hrp' & panel_outcome=='panel_positive' })
#   
#   fdr <- (function(){
#     tp <- l$hrd$abs[,'panel_positive']
#     fp <- l$hrp$abs[,'panel_positive']
#     fp / (tp+fp)
#   })()
# 
#   df_melt$bar_labels[fp_sel] <- with(df_melt[fp_sel,],{
#     labels <- paste0(
#       abs,
#       '\n(',signif(100*fdr,3),'%',')'
#     )
#     
#     labels[abs<15] <- paste0(labels[abs<15],'\n\n\n\n')
#     
#     return(labels)
#   })
#   
#   #df_melt$bar_labels[fp_sel] <- df_melt$abs[fp_sel]
#   
#   ## Force group ordering; set proper group names
#   group_labels <- c(
#     'hrp.panel_positive'='Patients incorrectly considered HRD',
#     'hrd.panel_negative'='HRD patients missed',
#     'hrd.panel_positive'='HRD patients detected'
#   )
#   
#   df_melt_split <- split(df_melt,df_melt$group)[group_order]
#   df_melt <- do.call(rbind, df_melt_split)
#   
#   df_melt <- as.data.frame(lapply(df_melt, function(i){
#     if(!is.numeric(i)){ i <- factor(i, unique(i)) }
#     return(i)
#   }))
# 
#   ggplot(df_melt, aes(x=panel_type, y=abs, fill=group, label=bar_labels)) +
#     geom_hline(yintercept=group_totals$hrd, size=0.25, linetype='dotted') +
#     geom_col(position=position_stack(), color='black', size=0.25) +
#     geom_text(position=position_stack(vjust=0.5), size=2.7) +
#     scale_fill_manual(values=c('#B1C4E4','#ede4d8','#F8D598'), labels=group_labels) +
#     
#     scale_y_continuous(breaks=seq(0, 50*round(group_totals$hrd/50), 50)) +
#     scale_x_discrete(expand=c(0.18, 0.25, 0.1, 0.25)) +
#     
#     labs(y='No. of patients', fill='Based on genetic screening:') +
#     #guides(fill=guide_legend(reverse=T, ncol=2)) +
#     
#     annotate('text', x=0.3, y=group_totals$hrd, size=3.25, angle=90, hjust=0, label='   (FDR)') +
#     annotate('text', x=0.3, y=group_totals$hrd/2, size=3.25, angle=90, hjust=0.5, label='(% of HRD patients)') +
#     
#     #theme_bw() +
#     theme_grey() +
#     theme(
#       panel.border=element_rect(fill=NA),
#       panel.background=element_rect(fill='#f1f1f1'),
#       panel.grid.major.x=element_blank(),
#       #panel.grid.major.y=element_blank(),
#       panel.grid.minor.y=element_blank(),
#       axis.title.x=element_blank()
#     )
# })()
# 
# pdf(paste0(base_dir,'/HMF_DR010_DR047/analysis/pancancer_overview/plots/panel_testing_confusion2.pdf'),9,4)
# grid.draw(plot_panel_confusion + theme(axis.text.x=element_blank()))
# dev.off()





# plotPanelConfusion <- function(
#   panel_names, col_values, legend_key_names, fill_colors=c('#B1C4E4','#F8D598')
# ){
#   ##
#   # panel_names=names(panel_outcomes$hrd)
#   # col_values=data.frame(
#   #   upper=length(hrd_samples) - panel_outcomes$hrd,
#   #   lower=panel_outcomes$hrd
#   # )
#   # legend_key_names=c(
#   #   'Is HRP based on genetic screening',
#   #   'Is HRD based on genetic screening'
#   # )
#   # fill_colors=c(
#   #   '#B1C4E4',
#   #   '#F8D598'
#   # )
#   
#   ## Sanity checks
#   if(length(panel_names) != nrow(col_values)){
#     stop('length(panel_names) != nrow(col_values)')
#   }
#   if(is.null(rownames(col_values))){
#     rownames(col_values) <- 1:nrow(col_values)
#   }
#   
#   ## Make bar labels
#   col_values_rel <- t(apply(as.data.frame(col_values),1,function(i){
#     i/sum(i)
#   }))
#   
#   bar_labels <- as.data.frame(do.call(rbind,lapply(1:nrow(col_values), function(i){
#     paste0(
#       signif(100*col_values_rel[i,],3),'%',
#       ' (',col_values[i,],')'
#     )
#   })))
#   colnames(bar_labels) <- names(col_values)
#   rownames(bar_labels) <- panel_names
#   
#   ## Prep data
#   plot_data <- list(
#     y=col_values_rel, 
#     labels=bar_labels
#   )
#   
#   plot_data2 <- lapply(plot_data, function(i){
#     i <- as.data.frame(i)
#     i <- cbind(panel_names=rownames(i),i)
#     rownames(i) <- NULL
#     melt(i, 'panel_names')
#   })
#   
#   plot_data2 <- Reduce(function(x,y){ merge(x,y,c('panel_names','variable'), sort=F) }, plot_data2)
#   colnames(plot_data2)[3:ncol(plot_data2)] <- names(plot_data)
#   
#   plot_data2 <- as.data.frame(lapply(plot_data2, function(i){
#     if(!is.numeric(i)){ i <- factor(i, unique(i)) }
#     return(i)
#   }))
#   
#   ggplot(plot_data2, aes(x=panel_names, y=y, fill=variable, label=labels)) +
#     geom_bar(stat='identity',position=position_stack(), color='black', size=0.25) +
#     geom_text(size=2.7, position=position_stack(vjust=0.5)) +
#     
#     scale_y_continuous(labels=function(x){ paste0(x*100,'%') }) +
#     scale_fill_manual(values=fill_colors, labels=legend_key_names) +
#     labs(y='% of patients', fill='Based on genetic screening:') +
#     
#     theme_grey() +
#     theme(
#       panel.border=element_rect(fill=NA, color='black'),
#       panel.grid.minor.y=element_blank(),
#       panel.grid.minor.x=element_blank(),
#       panel.grid.major.x=element_blank(),
#       
#       #legend.title=element_blank(),
#       #legend.spacing.x=unit(4,'pt'),
#       legend.key=element_rect(color=NA),
#       #legend.key.width=unit(1,'lines'),
#       #legend.key.height=unit(2,'lines'),
#       legend.justification=c(0,0.5),
#       
#       axis.text.x=element_text(angle=25, hjust=1, vjust=1),
#       axis.title.x=element_blank()
#     )
# }
# 
# ## Exec
# plots <- list()
# 
# plots$hrd_samples <- plotPanelConfusion(
#   panel_names=names(panel_outcomes$hrd),
#   
#   col_values=data.frame(
#     length(hrd_samples) - panel_outcomes$hrd,
#     panel_outcomes$hrd
#   ),
#   
#   legend_key_names=c(
#     'HRD patients missed',
#     'HRD patients detected'
#   ),
#   
#   fill_colors=c('#F8D598','#B1C4E4')
# )
# 
# 
# plots$fp <- plotPanelConfusion(
#   panel_names=names(panel_outcomes$hrd),
#   
#   col_values=data.frame(
#     panel_outcomes$hrp,
#     panel_outcomes$hrd
#   ),
#   
#   legend_key_names=c(
#     'Patients not correctly identified as HRD',
#     'Patients correctly identified as HRD'
#   ),
#   
#   fill_colors=c('#F7AEA9','#B1C4E4')
# )
# 
# plots <- lapply(plots, function(i){
#   i + theme(
#     axis.text.x=element_blank(),
#     axis.ticks.x=element_blank()
#   )
# })
# 
# ## Export
# plots_combined <- plot_grid(plotlist=plots, align='v', ncol=1)
# 
# pdf(paste0(base_dir,'/HMF_DR010_DR047/analysis/pancancer_overview/plots/panel_testing_confusion.pdf'),9,5)
# grid.draw(plots_combined)
# dev.off()



