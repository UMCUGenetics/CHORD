library(mutSigExtractor)

options(stringsAsFactors = F)

#========= Paths =========#
for(i in base_dir){
  if(dir.exists(i)){ 
    base_dir <- i 
    break
  }
}

#========= Determine MSI =========#
msi_samples <- (function(){
  contexts <- readRDS(paste0(base_dir,'/HMF_DR010_DR047/matrices/merged_contexts.rds'))
  contexts <- transformContexts(contexts, simplify.types=c('snv','indel','sv'))
  
  df <- as.data.frame(contexts[,c('del.rep','ins.rep')])
  df$indel.rep <- rowSums(df)
  
  df$has_msi <- df$indel.rep >= 14000
  
  rownames(df)[df$has_msi]
})()

#========= Main =========#
diplotypes_brca <- read.delim(paste0(base_dir,'/HMF_DR010_DR047/scripts/det_brca_status/hmf_gene_diplotypes_max_brca.txt.gz'))

detTrainingSet <- function(df, genes=c('BRCA2','BRCA1'), msi.samples){
  
  # df=diplotypes_brca
  # genes=c('BRCA2','BRCA1')
  # msi.samples=msi_samples

  ## Determine pre-conditions
  has_msi <- df$sample %in% msi.samples
  
  is_def <- df$hit_score >= 10 
  is_prof <- unlist(with(df,{
    Map(function(diplotype_origin, a1.max_score, a2.max_score){
      if(diplotype_origin=='cnv_germ'){
        a1.max_score==0 & a2.max_score <= 3

      } else if(diplotype_origin=='cnv_som'){
        a1.max_score==0 & a2.max_score <= 3

      } else if(diplotype_origin=='germ_som'){
        a1.max_score < 4 & a2.max_score < 3
      
      } else {
        FALSE
      }
    }, diplotype_origin, a1.max_score, a2.max_score, USE.NAMES=F)
  }))
  
  ## Determine in training set (def)
  in_training_set <- rep(FALSE,nrow(df))
  in_training_set[is_def & !has_msi] <- TRUE ## BRCA1/2 class
  #in_training_set[is_prof_confident] <- TRUE ## none class
  
  ## Combine info into one df
  out_pre <- cbind(
    df[c('sample','hgnc_symbol','hit_score','a1','a1.max_score','a2','a2.max_score')],
    has_msi,
    is_def,
    is_prof,
    in_training_set
  )
  
  #--------- Prof (part 1) ---------#
  l <- lapply(genes, function(i){ 
    out_pre[out_pre$hgnc_symbol==i,c('sample','is_prof')] 
  })
  names(l) <- genes
  
  df_det_is_prof_all <- Reduce(function(x,y){ merge(x,y,by='sample') },l)
  df_det_is_prof_all$is_prof_all <- apply(df_det_is_prof_all[,-1],1,all)
  
  #--------- Def + make final output ---------#
  ## Sort by genes order
  out <- do.call(rbind, lapply(genes, function(i){ 
    out_pre[out_pre$hgnc_symbol==i,] 
  }))
  
  ## Then, greedily resolve multiple hit cases
  out <- do.call(rbind, lapply(unique(out$sample),function(i){
    df_ss <- out[out$sample==i,]
    df_ss[which.max(df_ss$is_def),]
  }))
  
  ## Indicate which gene was deficient
  out$response <- with(out,{
    unlist(Map(function(is_def, hgnc_symbol){
      if(is_def){ hgnc_symbol }
      else{ 'none' }
    }, is_def, hgnc_symbol))
  })
  
  insertColAfter <- function(df, position, v, name){
    col_index <- which(position==colnames(df))
    df <- cbind(df[,1:col_index], v, df[,(col_index+1):ncol(df)])
    colnames(df)[col_index+1] <- name
    return(df)
  }
  
  #--------- Prof (part 2) ---------#
  out <- insertColAfter(
    out, 'is_prof',
    df_det_is_prof_all$is_prof_all[match(out$sample, df_det_is_prof_all$sample)],
    'is_prof_all'
  )
  
  out$in_training_set[out$is_prof_all] <- TRUE

  return(out)
}

annotation <- detTrainingSet(diplotypes_brca, genes=c('BRCA2','BRCA1'), msi_samples)

write.table(
  annotation, paste0(base_dir,'HMF_DR010_DR047/scripts/det_brca_status/hmf_brca_response.txt'),
  sep='\t',quote=F, row.names=F
)

