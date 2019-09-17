#' Predict the probability of homogolous recombination deficiency using mutational signatures
#' 
#' @description A wrapper for predict.randomForest() from the randomForest package
#'
#' @param sigs The output of extractSigsChord(), which is a dataframe containing the SNV, indel and
#' SV context counts.
#' @param rf.model The random forest model. Defaults to CHORD.
#' @param hrd.cutoff Samples greater or equal to this cutoff will be marked as HRD in the output
#' table of chordPredict(). Default is 0.5.
#' @param show.sigs Show the signature matrix in the output table?
#'
#' @return A dataframe containing per sample the probabilities of BRCA1/2 deficiency and HRD, as
#' well as the input signature values if show.sigs==T
#' @export
#'
#' @examples
#' ## Extract signatures
#' vcf_dir <- '/path_to_vcfs/'
#' vcf_snv <- paste0(vcf_dir,'SampleX_post_processed_v2.2.vcf.gz')
#' vcf_indel <- paste0(vcf_dir,'SampleX_post_processed_v2.2.vcf.gz')
#' vcf_sv <- paste0(vcf_dir,'SampleX_somaticSV_bpi.vcf.gz')
#' sigs <- extractSigsChord(vcf_snv, vcf_indel, vcf_sv, sample.name='SampleX')
#' 
#' ## Predict HRD probability with CHORD
#' chordPredict(sigs)

chordPredict <- function(sigs, rf.model=CHORD, hrd.cutoff=0.5, show.sigs=F){
  
  ## Converts the raw signature counts from extractSigsChord() to features used by CHORD
  sigs_split <- mutSigExtractor::splitDfRegex(sigs, c(snv='>',indel='[a-z]{3}[.]',sv='[A-Z]{3}'))
  sigs_processed <- mutSigExtractor::transformContexts(
    sigs_split,
    
    ## Simplify indels to types only (mh: flanking microhomology, rep: within repeat regions, none: 
    ## other indels) by ignoring indel length.
    ## Simplify SNVs to the 6 types of base substitions (C>A, C>G, C>T, T>A, T>C, T>G) by ignoring
    ## the immediate flanking nucleotides
    simplify.types = c('snv','indel'),
    
    ## Convert absolute counts to relative counts. Calculated per variant type.
    rel.types = c('snv','indel','sv')
  )
  
  df <- as.data.frame(predict(rf.model, sigs_processed, type='prob'))
  df <- df[,c('none','BRCA1','BRCA2')]
  colnames(df) <- paste0('p_',colnames(df))
  
  ## The probability of HRD is the probability of BRCA1 deficiency + BRCA2 deficiency
  df$p_hrd <- df$p_BRCA1 + df$p_BRCA2
  df$is_hrd <- df$p_hrd >= hrd.cutoff
  
  if(show.sigs){
    df <- merge(
      df,
      sigs_processed,
      by='row.names'
    )
  }
  
  df <- cbind(sample=rownames(df),df); rownames(df) <- NULL
  df <- df[order(df$p_hrd, decreasing=T),]
  
  return(df)
}



