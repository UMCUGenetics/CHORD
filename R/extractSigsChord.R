#' Extract signatures in the format compatible with CHORD
#' 
#' @description This function is a wrapper for the 3 functions from mutSigExtractor:
#' extractSigsSnv(), extractSigsIndel(), extractSigsSv(). Some post-processing is done to produce
#' compatible input for CHORD
#' 
#' @param vcf.snv Path to the vcf file containing SNVs
#' @param vcf.indel Path to the vcf file containing indels
#' @param vcf.sv Path to the vcf file containing SVs
#' @param sample.name The name of the sample as a character. Defaults to 'sample' if none is 
#' provided.
#' @param sv.caller SV vcfs are not standardized and therefore need to be parsed differently
#' depending on the caller. Currently supports 'manta' or 'gridss'.
#' @param output.path If a path is specified, the output is written to this path.
#' @param verbose Whether to print progress messages
#'
#' @return A 1-row data frame containing the mutational signature contributions
#' @export
#'
extractSigsChord <- function(
  vcf.snv, vcf.indel, vcf.sv, sample.name='sample', 
  sv.caller='gridss', output.path=NULL, verbose=T
){
  # vcf_dir='/Users/lnguyen//hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-DR047/data/160721_HMFregXXXXXXXX/'
  # vcf.snv=paste0(vcf_dir,'XXXXXXXX.vcf.gz')
  # vcf.indel=paste0(vcf_dir,'XXXXXXXX.vcf.gz')
  # vcf.sv=paste0(vcf_dir,'XXXXXXXX.purple.sv.ann.vcf.gz')
  
  sigs <- list()
  
  if(verbose){ message('Counting SNV trinucleotide contexts...') }
  sigs$snv <- extractSigsSnv(vcf.snv, vcf.filter='PASS', output='contexts', verbose=verbose)
  
  if(verbose){ message('Counting indel contexts (types x lengths)...') }
  sigs$indel <- extractSigsIndel(vcf.indel, vcf.filter='PASS', verbose=verbose)
  
  if(verbose){ message('Counting SV contexts (types x lengths)...') }
  sigs$sv <- extractSigsSv(
    vcf.sv, vcf.filter='PASS', sv.caller=sv.caller, output='contexts',
    sv.len.cutoffs = c(0, 10^3, 10^4, 10^5, 10^6, 10^7,Inf, verbose=verbose)
  )
  
  out <- do.call(cbind,lapply(sigs,t))
  rownames(out) <- sample.name
  
  if(is.null(output.path)){
    return(out)
  } else {
    write.table(out, output.path, sep='\t', quote=F)
  }
}


