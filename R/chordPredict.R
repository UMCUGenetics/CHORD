#' Predict the probability of homogolous recombination deficiency using mutational signatures
#' 
#' @description A wrapper for predict.randomForest() from the randomForest package
#'
#' @param features The output of extractSigsChord(), which is a dataframe containing the SNV, indel and
#' SV context counts.
#' @param rf.model The random forest model. Defaults to CHORD.
#' @param hrd.cutoff Samples greater or equal to this cutoff will be marked as HRD in the output
#' table of chordPredict(). Default is 0.5.
#' @param show.features Show the mutation context matrix in the output table?
#' @param verbose Show messages/warnings?
#'
#' @return A dataframe containing per sample the probabilities of BRCA1-type and BRCA2-type HRD, 
#' and HRD (= BRCA1-type HRD + BRCA2-type HRD)
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

chordPredict <- function(features, rf.model=CHORD, hrd.cutoff=0.5, show.features=F, verbose=T){
  
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
  
  ## Prediction
  df <- as.data.frame(predict(rf.model, features_processed, type='prob'))
  df <- df[,c('none','BRCA1','BRCA2')]
  colnames(df) <- paste0('p_',colnames(df))
  df$p_hrd <- df$p_BRCA1 + df$p_BRCA2
  
  ## QC
  df$has_msi <- with(features_split,{
    rowSums(indel[,grep('rep',colnames(indel)),drop=F]) > 14000
  })
  
  df$low_indel_load <- rowSums(features_split$indel) < 50
  df$low_sv_load <- rowSums(features_split$sv) < 30
  
  df$failed_qc <- with(df,{ has_msi | low_indel_load | low_sv_load })
  
  if(verbose & sum(df$failed_qc)>0){
    message(
      sum(df$failed_qc),' sample(s) failed QC:\n', 
      '  ', sum(df$has_msi),' with MSI (>14000 indels within repeats)\n',
      '  ', sum(df$low_indel_load), ' with <50 indels\n',
      '  ', sum(df$low_sv_load), ' with <30 SVs'
    )
  }
  
  ## Determine if sample is (confident) HRD
  df$is_hrd <- df$p_hrd >= hrd.cutoff
  df$is_hrd[ df$failed_qc ] <- NA
  
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



