#' Predict the probability of homogolous recombination deficiency using mutational signatures
#'
#' @description A wrapper for predict.randomForest() from the randomForest package
#'
#' @param features The output of extractSigsChord(), which is a dataframe containing the SNV, indel
#' and SV context counts.
#' @param rf.model The random forest model. Defaults to CHORD.
#' @param show.features Show the mutation context matrix in the output table?
#' @param hrd.cutoff Default=0.5. Samples greater or equal to this cutoff will be marked as HRD 
#' (is_hrd==TRUE).
#' @param min.indel.load Default=50. The minimum number of indels required to make an accurate HRD
#' prediction. Samples with fewer indels than this value will be marked as is_hrd==NA (HR status 
#' could not be confidently determined).
#' @param min.sv.load Default=30. The minimum number of SVs required to make an accurate prediction
#' of BRCA1-type vs. BRCA2-type HRD. Samples with fewer SVs than this value will be marked as 
#' hrd_type==NA (HRD type could not be confidently determined).
#' @param min.msi.indel.rep Default=14000 (changing this value is not advised). Samples with more 
#' indels within repeats than this threshold will be considered to have microsatellite instability.
#' @param detailed.remarks If TRUE, shows min.indel.load and min.sv.load numbers in the remarks columns
#' @param verbose Show messages/warnings?
#'
#' @return A dataframe containing per sample the probabilities of BRCA1-type and BRCA2-type HRD, and
#' HRD (= BRCA1-type HRD + BRCA2-type HRD)
#' @export
#' 
#' @examples
#' ## Extract mutation contexts
#' vcf_dir <- '/path_to_vcfs/'
#' vcf_snv <- paste0(vcf_dir,'SampleX_post_processed_v2.2.vcf.gz')
#' vcf_indel <- paste0(vcf_dir,'SampleX_post_processed_v2.2.vcf.gz')
#' vcf_sv <- paste0(vcf_dir,'SampleX_somaticSV_bpi.vcf.gz')
#' contexts <- extractSigsChord(vcf_snv, vcf_indel, vcf_sv, sample.name='SampleX')
#' 
#' ## Predict HRD probability with CHORD
#' chordPredict(contexts)

chordPredict <- function(
  features, rf.model=CHORD, show.features=F,
  hrd.cutoff=0.5, min.indel.load=50, min.sv.load=30, min.msi.indel.rep=14000,
  detailed.remarks=F, verbose=T
){
  
  ## Converts the raw signature counts from extractSigsChord() to features used by CHORD
  features_split <- mutSigExtractor::splitDfRegex(features, c(snv='>',indel='[a-z]{3}[.]',sv='[A-Z]{3}'))
  features_processed <- mutSigExtractor::transformContexts(
    features_split,
    
    ## Simplify indels to types only (mh: flanking microhomology, rep: within repeat regions, none: 
    ## other indels) by ignoring indel length.
    ## Simplify SNVs to the 6 types of base substitions (C>A, C>G, C>T, T>A, T>C, T>G) by ignoring
    ## the immediate flanking nucleotides
    simplify.types = c('snv','indel'),
    
    ## Convert absolute counts to relative counts. Calculated per variant type.
    rel.types = c('snv','indel','sv')
  )
  
  #--------- Prediction ---------#
  df <- as.data.frame(predict(rf.model, features_processed, type='prob'))
  df <- df[,c('none','BRCA1','BRCA2')]
  colnames(df) <- paste0('p_',colnames(df))
  df$p_hrd <- df$p_BRCA1 + df$p_BRCA2
  
  #--------- QC ---------#
  qc <- list()
  qc$has_msi <- with(features_split,{
    rowSums(indel[,grep('rep',colnames(indel)),drop=F]) > min.msi.indel.rep
  })
  
  qc$low_indel_load <- rowSums(features_split$indel) < min.indel.load
  qc$low_sv_load <- rowSums(features_split$sv) < min.sv.load
  qc <- as.data.frame(qc)
  
  failed_qc <- with(qc,{ has_msi | low_indel_load | low_sv_load })
  
  if(verbose & sum(failed_qc)>0){
    
    qc_messages <- list(
      has_msi=paste0('  Critical: ', sum(qc$has_msi),' with MSI (>',min.msi.indel.rep,' indels within repeats)\n'),
      low_indel_load=paste0('  Critical: ',sum(qc$low_indel_load),' with <', min.indel.load,' indels\n'),
      low_sv_load=paste0('  Non-critical: ', sum(qc$low_sv_load), ' with <', min.sv.load, ' SVs')
    )
    
    message(
      sum(failed_qc),' sample(s) failed QC:\n', 
      if(sum(qc$has_msi)>0){ qc_messages$has_msi } else { '' },
      if(sum(qc$low_indel_load)>0){ qc_messages$low_indel_load } else { '' },
      if(sum(qc$low_sv_load)){ qc_messages$low_sv_load } else { '' }
    )
  }
  
  ## Informative QC tags
  if(!detailed.remarks){
    failed_qc_strings <- colnames(qc)
    names(failed_qc_strings) <- colnames(qc)
  } else {
    ## Note to self: make sure order of qc names is correct upon editing
    failed_qc_strings <- c(
      has_msi=paste0('Has MSI (>',min.msi.indel.rep,' indel.rep)'),
      low_indel_load=paste0('<',min.indel.load,' indels'),
      low_sv_load=paste0('<',min.sv.load,' SVs')
    )
  }
  
  ## is_hrd qc tags
  df_qc_is_hrd <- data.frame(
    has_msi=ifelse(qc$has_msi,failed_qc_strings['has_msi'],''),
    low_indel_load=ifelse(qc$low_indel_load,failed_qc_strings['low_indel_load'],'')
  )
  
  qc_is_hrd <- with(df_qc_is_hrd,{
    paste(has_msi, low_indel_load,sep=';')
  })
  qc_is_hrd <- gsub('^;|;$','',qc_is_hrd)
  qc_is_hrd[nchar(qc_is_hrd)==0] <- ''
  
  ## hrd_type qc tags
  qc_hrd_type <- ifelse(
    qc$low_sv_load,
    failed_qc_strings['low_sv_load'],
    ''
  )
  
  qc_out <- data.frame(
    remarks_hr_status=qc_is_hrd, 
    remarks_hrd_type=qc_hrd_type
  )
  
  #--------- Determine if sample is HRD (only if sample has enough indels) ---------#
  df$hr_status <- ifelse(
    df$p_hrd >= hrd.cutoff,
    'HR_deficient','HR_proficient'
  )
  df$hr_status[ qc$low_indel_load | qc$has_msi ] <- 'cannot_be_determined'
  
  #--------- Determine HRD type ---------#
  df$hrd_type <- unlist(
    Map(function(p_BRCA1, p_BRCA2, p_hrd){
      if(p_hrd>=hrd.cutoff){
        c('BRCA1_type','BRCA2_type')[ which.max(c(p_BRCA1,p_BRCA2)) ]
      } else {
        'none'
      }
    }, df$p_BRCA1, df$p_BRCA2, df$p_hrd)
  )
  df$hrd_type[ qc$low_sv_load | qc$low_indel_load ] <- 'cannot_be_determined'
  
  df <- cbind(df, qc_out)
  
  if(show.features){
    df <- merge(
      cbind(sample=rownames(df),df),
      cbind(sample=rownames(df),features_processed),
      by='sample'
    )
  } else {
    df <- cbind(sample=rownames(df),df); rownames(df) <- NULL
  }
  
  return(df)
}



