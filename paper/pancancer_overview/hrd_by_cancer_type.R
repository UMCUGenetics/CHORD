library(RColorBrewer)
library(reshape2)
library(grid)
library(gridExtra)
library(ggplot2)
library(cowplot)
library(magrittr)

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
rank_order_clust <- readRDS(paste0(base_dir,'/HMF_DR010_DR047/analysis/pancancer_overview/data/rank_order_clust.rds'))
l_m_diplotypes <- readRDS(paste0(base_dir,'/HMF_DR010_DR047/analysis/pancancer_overview/data/l_m_diplotypes.rds'))

## Rm duplicate biopsies and keep samples with >=50 indels
analysis_sample_selection <- read.delim(paste0(base_dir,'/HMF_DR010_DR047/analysis/pancancer_overview/data/analysis_sample_selection.txt'))
selected_samples <- analysis_sample_selection[analysis_sample_selection$is_selected,'sample_id']

#========= Master dataframe =========#
summary_hrd <- (function(){
   
   #--------- Tag samples with genetic cause ---------#
   df <- data.frame(
      sample=names(rank_order_clust$clusters),
      cluster=rank_order_clust$clusters,
      is_hrd=TRUE,
      row.names=NULL
   )
   df$genotype <- names(rank_order_clust$cluster_names)[ match(df$cluster,rank_order_clust$cluster_names) ]
   
   df$genotype_monoall <- (function(){
      m <- rank_order_clust$df_ranked
      gene_cols <- grep('^p_',colnames(m),invert=T)
      apply(m[,gene_cols], 1, function(i){
         colnames(m[,gene_cols])[ which.max(i) ]
      })
   })()
   
   #--------- Add biallelic effect ---------#
   flattenDiplotypeMatrix <- function(m){
      #m=l_m_diplotypes$eff
      
      genes <- unique(sapply(strsplit(colnames(m),'_'),`[`,1))
      gene_cols <- lapply(genes, function(i){ 
         grep(i, colnames(m))
      }); names(gene_cols) <- genes
      
      ## Reorder eff matrix to match with df
      if(!all(df$sample %in% rownames(m))){
         warning("df$sample and rownames(m) are not completely intersecting")
      }
      
      m <- m[match(df$sample,rownames(m)),]
      
      ## Select diplotype (value pairs) for correct gene
      m$genotype <- df$genotype_monoall
      # m$genotype <- gsub('unknown_', '', m$genotype)
      # m$genotype <- gsub('-type', '', m$genotype)
      
      m_flat <- t(apply(as.matrix(m),1,function(i){
         genotype <- i['genotype']
         if(genotype %in% genes){
            i[ gene_cols[[genotype]] ]
         } else {
            c('none','none')
         }
      })); colnames(m_flat) <- c('a1','a2')
      
      return(m_flat)
   }
   
   m_diplotypes_flat <- (function(){
      m_eff <- flattenDiplotypeMatrix(l_m_diplotypes$eff)
      
      m_a_origin <- flattenDiplotypeMatrix(l_m_diplotypes$a_origin)
      colnames(m_a_origin) <- paste0(colnames(m_a_origin),'.origin')
      
      m_a_score <- flattenDiplotypeMatrix(l_m_diplotypes$score)
      colnames(m_a_score) <- paste0(colnames(m_a_score),'.score')
      
      data.frame(m_eff, m_a_origin, m_a_score)
   })()
   
   df$diplotype_origin <- apply(m_diplotypes_flat,1,function(i){
      a1 <- i[1]
      a2 <- i[2]
      a1.origin <- i[3]
      a2.origin <- i[4]
      
      if(a1 %in% c('full_gene_loss','trunc')){
         a1
      } else if(a1=='loh'){
         if(a2!='none'){
            paste0('loh+',a2.origin)
         } else {
            'loh+unknown'
         }
      } else if(a1!='none' & a2!='none'){
         paste0(a1.origin,'+',a2.origin)
      } else {
         'unknown'
      }
   })
   
   ## Format diplotype origin
   df <- within(df,{
      diplotype_origin <- gsub('[+]', ' + ', diplotype_origin)
      diplotype_origin <- gsub('full_gene_loss', 'Deep deletion', diplotype_origin)
      diplotype_origin <- gsub('trunc', 'Deep deletion', diplotype_origin)
      diplotype_origin <- gsub('loh', 'LOH', diplotype_origin)
   })
   
   ## Format genotypes (pt.2)
   df$genotype <- gsub('unknown_BRCA1-type','unkn. (BRCA1-type HRD)',df$genotype)
   df$genotype <- gsub('unknown_BRCA2-type','unkn. (BRCA2-type HRD)',df$genotype)
   
   #--------- Add HRP samples to df ---------#
   df_prof <- df[1,]
   counter <- 0
   for(i in df){
      counter <- counter+1
      if(is.numeric(i)){ df_prof[counter] <- 0 }
      else if(is.logical(i)){ df_prof[counter] <- FALSE }
      else { df_prof[counter] <- 'none' }
   }
   
   prof_samples <- selected_samples[ !(selected_samples %in% df$sample) ]
   df_prof <- df_prof[rep(1,length(prof_samples)),]
   df_prof$sample <- prof_samples
   
   df <- rbind(df, df_prof)
   
   if(!all(df$sample %in% analysis_sample_selection$sample_id)){
      warning("df$sample and metadata$sample_id are not completely intersecting")
   }
   
   #--------- Add cancer type ---------#
   df$cancer_type <- analysis_sample_selection$primary_tumor_location[ match(df$sample,analysis_sample_selection$sample_id) ]
   
   return(df)
})()


####################################################################################################
# HRD cause by cancer type                                                                         #
####################################################################################################

#========= Main counting function =========#
countFeature1ByFeature2 <- function(
   df, feature.split, feature.count, rm.rownames=T, 
   calc.rel.counts=NULL, 
   calc.unsplit.counts=F, unsplit.name='Pan-cancer'
){
   # df=summary_hrd
   # feature.split='cancer_type'
   # feature.count='is_hrd'
   # row.order=hrd_by_cancer_type$feature_split
   
   df_split <- split(df, df[,feature.split])
   df_split <- df_split[unique(df[,feature.split])] ## preserve original order
   
   ## Main
   values_uniq <- unique(df[,feature.count])
   values_uniq <- structure(rep(0,length(values_uniq)),names=values_uniq)
   out <- do.call(rbind, lapply(df_split, function(i){
      #i=df_split[[1]]
      tab <- table(i[,feature.count])
      v <- values_uniq
      v[names(tab)] <- tab
      return(v)
   }))
   
   out <- as.data.frame(out)
   
   ## Calculate pancancer counts
   if(calc.unsplit.counts){
      tab <- table(df[,feature.count])
      v <- values_uniq
      v[names(tab)] <- tab
      out <- rbind(v, out)
      rownames(out)[1] <- 
         if(is.null(unsplit.name)){ 'All' }
         else { unsplit.name }
   }
   
   if(!is.null(calc.rel.counts)){
      if(calc.rel.counts=='macro'){ out <- out/nrow(df) }
      if(calc.rel.counts=='micro'){ out <- out/rowSums(out) }
   }
   
   ## feature.split as column
   if(rm.rownames){
      out <- cbind(feature_split=rownames(out),out)
      rownames(out) <- NULL
   }
   
   return(out)
}


#========= HRD by cancer type =========#
calcHrdByCancerType <- function(summary.hrd, min.hrd.freq=5){
   #summary.hrd=summary_hrd
   #summary.hrd=subset(summary_hrd, cluster %in% c(1,2,3,5))
   
   df <- as.data.frame(
      countFeature1ByFeature2(summary.hrd, 'cancer_type', 'is_hrd', rm.rownames=F, calc.unsplit.counts=T)
   )
   df$total <- df[['FALSE']] + df[['TRUE']]
   colnames(df)[1:2] <- c('abs','rel')
   df$rel <- df$abs/df$total
   
   df <- data.frame(feature_split=rownames(df), df)
   
   ## 
   df$split_var <- ifelse(df$abs <= min.hrd.freq, 'low', 'high')
   df$split_var[df$feature_split=='Pan-cancer'] <- 'except'
   
   ## Deal split sort cancer types with low/high absolute number of hrd samples
   df_split <- split(df, df$split_var)
   out <- do.call(rbind, lapply(df_split, function(i){
      i[order(i$rel, decreasing=T),]
   }))
   out$feature_split <- factor(out$feature_split,out$feature_split)
   out$split_var <- NULL
   
   ## Store metadata
   class(out) <- c(class(out), min.hrd.freq=min.hrd.freq)
   rownames(out) <- NULL
   
   return(out)
}

hrd_by_cancer_type <- calcHrdByCancerType(summary_hrd)

plotHrdByCancerType <- function(
   df, show.min.hrd.freq.line=T, min.hrd.freq.line.hjust=0,
   feature.name=NULL, flip.category.axis=T,
   y.upper.lim.scale=1.5
){
   #df=hrd_by_cancer_type
   
   min_hrd_freq <- as.integer(class(df)['min.hrd.freq'])
   df$feature_split <- factor(df$feature_split,df$feature_split)
   
   df$label <- with(df,{
      perc <- paste0( round(100*rel,1), '%' )
      #label <- perc
      #label <- paste0(perc, ' (', abs,')')
      #label <- paste0(perc, ' (', abs,'/',total,')')
      label <- paste0(abs,' / ',total)
      label[abs==0 & total==0] <- ''
      return(label)
   })
   
   bar_fills <- rep('#51A152', length(unique(df$feature_split)))
   #bar_fills[df$abs < min_hrd_freq] <- '#bcdcbc'
   bar_fills[1] <- '#C82A26'
   
   bar_colors <- rep('black', length(unique(df$feature_split)))
   bar_colors[df$abs==0] <- NA
   
   plot <- ggplot(df, aes(x=feature_split, y=rel)) + 
      geom_bar(stat='identity', size=0.25, width=0.85, fill=bar_fills, color=bar_colors) +
      geom_text(aes(label=label), angle=90, hjust=-0.1, size=2.5) +
      #geom_text(aes(label=label,y=-0.08), angle=90, hjust=0, size=2.5) +
      scale_y_continuous(
         labels=scales::percent_format(accuracy=1),
         #limits=c(-0.08, NA)
         limits=c(0, max(df$rel)*y.upper.lim.scale)
      ) +
      theme_bw() +
      theme(
         #panel.grid.major.x=element_blank(),
         #panel.grid.minor.x=element_blank(),
         #panel.grid.minor.y=element_blank(),
         axis.line=element_blank(),
         axis.title.x=element_blank(),
         axis.text.x.bottom=element_text(angle=90,hjust=1,vjust=0.5),
         axis.text.x.top=element_text(angle=90,hjust=0,vjust=0.5),
         #axis.ticks.x=element_blank(),
         axis.title.y=element_text(size=9)
      )
   
   if(show.min.hrd.freq.line){
      ann_line_x <- which.max(df$abs <= min_hrd_freq) + min.hrd.freq.line.hjust
      plot <- plot +
         geom_vline(xintercept=ann_line_x, linetype='dashed', size=0.2) +
         annotate(
            'text', x=ann_line_x, y=max(df$rel*y.upper.lim.scale), 
            #label=sprintf("HRD freq. < %s",min_hrd_freq),
            label=paste0("'HRD freq.' <= ",min_hrd_freq), parse=T,
            hjust=-0.1, vjust=2, size=2.5
         )
   }
   
   if(flip.category.axis){ plot <- plot + scale_x_discrete(position='top') }
   if(!is.null(feature.name)){ plot <- plot + ylab(feature.name) }
   
   return(plot)
}


#========= Count features by cancer type =========#
#--------- Counts ---------#
sel_features <- c('genotype','diplotype_origin')
feature_counts <- lapply(sel_features, function(i){
   out <- countFeature1ByFeature2(
      summary_hrd,
      feature.split='cancer_type',
      feature.count=i,
      calc.unsplit.counts=T
   )
   out <- out[ match(hrd_by_cancer_type$feature_split, out$feature_split), ]
   out$feature_split <- factor(out$feature_split,out$feature_split)
   out <- out[colnames(out)!='none']
})
names(feature_counts) <- sel_features

#--------- Simplify LOH + mut groups  ---------#
feature_counts$biall_hit_type <- within(feature_counts$diplotype_origin,{
   #`LOH + small mut.` <- `LOH + germ` + `LOH + som`
   # `LOH + germline mut.` <- `LOH + germ`
   # `LOH + somatic mut.` <- `LOH + som`
   # `LOH + germ` <- `LOH + som` <- NULL
   
   `2x small mut.` <- `germ + som`
   `germ + som` <- NULL
})

colnames(feature_counts$biall_hit_type)[3:4] <- c('LOH + somatic mut.','LOH + germline mut.')


#--------- Biallelic hit origin  ---------#
feature_counts$mut_origin <- with(feature_counts$diplotype_origin,{
  `Germline + somatic` <- `LOH + germ` + `germ + som`
  `2x somatic` <- `Deep deletion` + `LOH + som`
  Unknown <- `LOH + unknown` + unknown
   
  data.frame(feature_split, `Germline + somatic`, `2x somatic`, Unknown, check.names=F)
})

#========= Export counts =========#
calcRelCounts <- function(df){
   #df=feature_counts$biall_hit_type
   rel_counts <- df[,-1]/rowSums(df[,-1])
   cbind(feature_split=df[,1], rel_counts)
}

library(openxlsx)

l_xlsx <- list()
l_xlsx[['Frequency of HRD']] <- hrd_by_cancer_type

l_xlsx[['Gene deficiency']] <- feature_counts$genotype
l_xlsx[['Gene deficiency (rel.)']] <- calcRelCounts(feature_counts$genotype)

l_xlsx[['Biall. hit type']] <- feature_counts$biall_hit_type
l_xlsx[['Biall. hit type (rel.)']] <- calcRelCounts(feature_counts$biall_hit_type)

l_xlsx[['Origin of biall. hit']] <- feature_counts$mut_origin
l_xlsx[['Origin of biall. hit (rel.)']] <- calcRelCounts(feature_counts$mut_origin)

formatXlsxTable <- function(df){
   # isNotFloatingPoint <- function(x){
   #    all(!grepl('^\\d+[.]',x))
   # }
   
   out <- as.data.frame(lapply(df,function(i){
      if(is.character(i)){ i }
      else if(is.numeric(i)){
         i[is.na(i)] <- 0
         i
      }
      else { i }
   }), check.names=F)
   
   colnames(out)[1] <- 'Cancer type'
   return(out)
}

l_xlsx <- lapply(l_xlsx, formatXlsxTable)

write.xlsx(
   l_xlsx, paste0(base_dir,'/HMF_DR010_DR047/analysis/pancancer_overview/supp_tables/hrd_by_cancer_type.xlsx')
)


#--------- Plot parameters ---------#
## Order of elements will be enforced in the plots
group_fills <- list()

group_fills$genotype <- list(
   BRCA1_type=c(
      'BRCA1'='#2C867C', ## Blue
      'unkn. (BRCA1-type HRD)'='#BCE6DF' ## Light blue
   ),

   BRCA2_type=c(
      'BRCA2'='#6c3809',#'#783f0b' ## Dark brown
      'RAD51C'='#B06E23',#'#b04b23', ## brown
      'PALB2'='#D7B66A', ## light brown
      'unkn. (BRCA2-type HRD)'='#f7eed3' ## pale brown
   )
)

group_fills$biall_hit_type <- c(
   'unknown'='white',
   '2x small mut.'='#FEF7C2', ## light yellow

   'LOH + unknown'='#d1c0df',
   #'LOH + small mut.'='#663096', ## purple
   'LOH + germline mut.'='#8459ab', ## light purple
   'LOH + somatic mut.'='#512678', ## purple

   'Deep deletion'='#ff42a8' ##pink
)

group_fills$mut_origin <- c(
   Unknown='white',
   `2x somatic`='#8DB7D3',
   `Germline + somatic`='#3163A0'
)

#========= Plotting =========#
plotFeatureByFeature <- function(
   df, group.fills=NULL, category.alphas=rep(1, nrow(df)),
   feature.name=NULL, 
   legend.position=c(0.7,0.5), legend.justification=c(0,0.5),
   hide.category.axis=F, geom.bar.position='fill'
){
   # df=feature_counts$mut_origin
   # group.fills=group_fills$mut_origin
   # feature.name='Mut. origin'
   
   df$feature_split <- factor(df$feature_split,df$feature_split)
   
   df_melt <- melt(df,'feature_split')
   colnames(df_melt) <- c('feature_split','feature','count')
   
   ## Force ordering and assign colors based on group.fills
   if(!is.null(group.fills)){
      df_melt$feature <- factor(df_melt$feature,names(group.fills))
   }

   ## Main
   plot <- ggplot(df_melt, aes(feature_split, count, group=feature, fill=feature, alpha=feature_split)) +
      geom_bar(stat='identity', position=geom.bar.position,color='black', width=0.8, size=0.25) +
      scale_fill_manual(values=group.fills, labels=function(x){ substr(x, 1, 1) <- toupper(substr(x, 1, 1)); x }) +
      scale_alpha_manual(values=category.alphas, guide=F) +
      scale_y_continuous(labels=scales::percent_format(accuracy=1)) +
      
      theme_bw() +
      
      theme(
         panel.border=element_rect(fill=NA),
         panel.grid.minor=element_blank(),
         axis.line=element_blank(),
         axis.title.x=element_blank(),
         axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
         axis.title.y=element_text(size=9),
         #legend.title=element_blank(), ## For some reason, this still leaves some white space
         legend.title=element_text(size=0),
         legend.spacing.x=unit(0.1, 'cm'),
         legend.key.size=unit(0.6,"line"),
         legend.text=element_text(size=7),
         legend.box.background=element_rect(fill=NA)
      )
   
   if(!is.null(feature.name)){ plot <- plot + ylab(feature.name) }
   if(!is.null(legend.position)){ plot <- plot + theme(legend.position=legend.position) }
   if(!is.null(legend.position)){ plot <- plot + theme(legend.justification=legend.justification) }
   
   if(hide.category.axis){
      plot <- plot + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
   }
   
   return(plot)
}

#========= Combine =========#
#--------- Add dummy row to insert space between high/low abs HRD incidence cancer types ---------#
insertRow <- function(df, df.ins, after.row){
   df_top <- df[1:after.row,]
   df_bottom <- df[(after.row+1):nrow(df),]
   
   colnames(df.ins) <- colnames(df)
   rbind(df_top, df.ins, df_bottom)
}

addDummyRows <- function(df, after.rows){
   #df=hrd_by_cancer_type
   #df=feature_counts$genotype
   
   counter <- 0
   
   for(i in sort(after.rows, decreasing=T)){
      #i=1
      counter <- counter + 1
      row_insert <- as.data.frame(lapply(df, function(j){
         if(is.character(j) || is.factor(j)){
            paste(rep(' ', counter), collapse='')
         } else {
            0
         }
      }), )
      colnames(row_insert) <- colnames(df)
      df <- insertRow(df, row_insert, i)
   }
   
   return(df)
}

row_insert_positions <- c(7)
custom_data <- list(
   incidence=addDummyRows(hrd_by_cancer_type, row_insert_positions),
   genotype=addDummyRows(feature_counts$genotype, row_insert_positions),
   biall_hit_type=addDummyRows(feature_counts$biall_hit_type, row_insert_positions),
   mut_origin=addDummyRows(feature_counts$mut_origin, row_insert_positions)
)

#--------- Custom genotype plot ---------#
genotype_plot <- (function(){
   
   ## Hacky solution to split up BRCA2/BRCA1-type HRD keys and add title
   fill_pal <- with(group_fills$genotype, c(BRCA1_type, rev(BRCA2_type)))
   legend_breaks <- with(group_fills$genotype, c(names(BRCA2_type), names(BRCA1_type)))
   
   legend_labels <- legend_breaks
   legend_labels[grep('^unkn',legend_labels)] <- 'Unknown             '
   
   plot <- plotFeatureByFeature(
      custom_data$genotype,
      group.fills = fill_pal,
      feature.name='Gene deficiency',
      hide.category.axis=T
   )
   
   plot <- plot + 
      scale_fill_manual(
         breaks=legend_breaks,
         labels=legend_labels,
         values=fill_pal
      ) +
      guides(fill=guide_legend(title='BRCA2-type HRD    BRCA1-type HRD',ncol=2, nrow=4)) +
      theme(legend.title=element_text(hjust=0, size=7, face='bold'))
   
   return(plot)
})()

#--------- Main ---------#
plots <- list(
   incidence=plotHrdByCancerType(
      custom_data$incidence, 
      feature.name='Freq. HRD', y.upper.lim.scale=1.9
   ),
   
   genotype=genotype_plot + geom_vline(xintercept=8, linetype='dashed', size=0.2),
   
   biall_hit_type=plotFeatureByFeature(
      custom_data$biall_hit_type,
      group.fills = group_fills$biall_hit_type,
      feature.name='Biallelic hit type',
      hide.category.axis=T
   ) + geom_vline(xintercept=8, linetype='dashed', size=0.2),
   
   mut_origin=plotFeatureByFeature(
      custom_data$mut_origin,
      group.fills = group_fills$mut_origin,
      feature.name='Origin of biallelic hit',
      hide.category.axis=F
   ) + geom_vline(xintercept=8, linetype='dashed', size=0.2)
)

#--------- Format x axis ---------#
## Set dummy cancer type tick to white for first/last plot
tick_colors <- rep('black',nrow(custom_data$incidence))
tick_colors[custom_data$incidence$feature_split==' '] <- 'white'

## Make 'Pancancer' bold
axis_text_x_face <- rep('plain',nrow(custom_data$incidence))
axis_text_x_face[1] <- 'bold'

x_axis_theme <- theme(
   axis.ticks.x=element_line(color=tick_colors), 
   axis.text.x=element_text(face=axis_text_x_face)
)

plots[[1]] <- plots[[1]] + x_axis_theme
plots[[length(plots)]] <- plots[[length(plots)]] + x_axis_theme

plots_combined <- plot_grid(
   plotlist=plots, align='v', axis='tblr', ncol=1,
   rel_heights=c(1.75, 1, 1, 1.75)
)

pdf(paste0(base_dir,'/HMF_DR010_DR047/analysis/pancancer_overview/plots/hrd_cause_by_cancer_type.pdf'),8,7)
grid.draw(plots_combined)
dev.off()


####################################################################################################
# Cancer types, biallelic loss clusters vs. non-biallelic loss clusters                            #
####################################################################################################
cancer_type_order <- hrd_by_cancer_type$feature_split[hrd_by_cancer_type$feature_split %in% summary_hrd$cancer_type]

cancer_type_counts_split <- (function(){
   l <- list(
      `Biallelic loss\n(clusters: 1, 2, 3, 5)`=table(subset(summary_hrd, cluster %in% c(1,2,3,5))$cancer_type),
      `Non-biallelic loss\n(clusters: 4, 6)`=table(subset(summary_hrd, cluster %in% c(4,6))$cancer_type)
   )
   
   l <- lapply(l, function(i){
      #i <- l[[1]]
      missing_cancer_types <- unique(summary_hrd$cancer_type)[ !(unique(summary_hrd$cancer_type) %in% names(i)) ]
      i[missing_cancer_types] <- 0
      as.table(i)
   })
   
   pd <- do.call(rbind, lapply(names(l), function(i){
      #i <- names(l)[[1]]
      df <- data.frame(
         cancer_type=names(l[[i]]),
         abs=as.vector(l[[i]])
      )
      
      df <- df[match(cancer_type_order, df$cancer_type),]
      
      df$cancer_type <- as.character(df$cancer_type)
      df$cancer_type[!(df$cancer_type %in% c('Ovary','Pancreas','Prostate','Breast','Biliary','Urinary tract'))] <- 'Other'
      
      df_split <- split(df, df$cancer_type=='Other')
      
      df_split[['TRUE']] <- data.frame(
         cancer_type='Other',
         abs=sum(df_split[['TRUE']]$abs)
      )
      
      df <- do.call(rbind, df_split)
      
      df$rel <- df$abs/sum(df$abs)
      df$label <- paste0(
         df$abs,'\n',
         '(',signif(df$rel*100,2),'%)'
      )
      
      df$group <- i
      
      
      return(df)
   }))
   
   as.data.frame(lapply(pd, function(i){
      if(is.numeric(i)){ i }
      else { factor(i, unique(i)) }
   }))
})()

cancer_type_counts_split_plot <- 
   ggplot(cancer_type_counts_split, aes(cancer_type, abs, group=group, label=label)) +
      geom_bar(stat='identity') +
      geom_text(vjust=-0.5, size=2.5) +
      facet_grid(group~.) +
      ylab('Number of patients') +
      ylim(0,47) +
      theme_bw() +
      theme(
         axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
         axis.title.x=element_blank()
      )

pdf(paste0(base_dir,'/HMF_DR010_DR047/analysis/pancancer_overview/plots/cancer_type_counts_yes_vs_no_biall_loss.pdf'),7,5)
grid.draw(cancer_type_counts_split_plot)
dev.off()
