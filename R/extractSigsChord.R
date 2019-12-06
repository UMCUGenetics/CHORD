#' Extract mutation contexts in the format compatible with CHORD
#' 
#' @description This function is a wrapper for the 3 functions from mutSigExtractor:
#' extractSigsSnv(), extractSigsIndel(), extractSigsSv(). Some post-processing is done to produce
#' compatible input for CHORD
#' 
#' @param vcf.snv Path to the vcf file containing SNVs
#' @param vcf.indel Path to the vcf file containing indels. By default vcf.indel==vcf.snv
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
  vcf.snv, 
  vcf.indel=vcf.snv, 
  vcf.sv,
  sample.name='sample',
  vcf.filters=c(snv=NA,indel=NA,sv=NA),
  sv.caller='gridss', output.path=NULL, verbose=F,
  ref.genome=DEFAULT_GENOME
){
  
  ######### Load vcfs #########
  if(verbose){ message('\n#====== Loading variants from vcfs ======#') }
  
  variants <- list()
  
  if(verbose){ message('\n## SNVs') }
  variants$snv <- variantsFromVcf(
    vcf.snv, mode='snv_indel', 
    vcf.filter=vcf.filters['snv'],
    ref.genome=ref.genome,
    verbose=verbose
  )
  
  if(verbose){ message('\n## Indels') }
  if(vcf.indel==vcf.snv){
    if(verbose){ message('vcf file is the same for both SNVs and indels. Skipping reading vcf for indels') }
    variants$indel <- variants$snv
  } else {
    variants$indel <- variantsFromVcf(
      vcf.indel, mode='snv_indel', 
      vcf.filter=vcf.filters['indel'], 
      ref.genome=ref.genome,
      verbose=verbose
    )
  }
  
  if(verbose){ message('\n## SVs') }
  variants$sv <- variantsFromVcf(
    vcf.sv, mode='sv', 
    vcf.filter=vcf.filters['sv'], 
    sv.caller=sv.caller
  )
  
  ######### Count contexts #########
  if(verbose){ message('\n#====== Counting mutation contexts ======#') }
  sigs <- list()
  
  if(verbose){ message('\n## Single base substitutions') }
  sigs$snv <- extractSigsSnv(bed=variants$snv, output='contexts', ref.genome=ref.genome, verbose=verbose)
  
  if(verbose){ message('\n## Indel contexts (types x lengths)') }
  sigs$indel <- extractSigsIndel(bed=variants$indel, ref.genome=ref.genome, verbose=verbose)
  
  if(verbose){ message('\n## SV contexts (types x lengths)') }
  sigs$sv <- extractSigsSv(
    df=variants$sv, sv.caller=sv.caller, output='contexts',
    sv.len.cutoffs = c(0, 10^3, 10^4, 10^5, 10^6, 10^7,Inf, verbose=verbose)
  )
  
  ######### Export #########
  if(verbose){ message('\n#====== Exporting output =========#') }
  out <- do.call(cbind,lapply(sigs,t))
  rownames(out) <- sample.name
  
  if(is.null(output.path)){
    if(verbose){ message('output.path not specified. Directly returning output') }
    return(out)
  } else {
    if(verbose){ message('Writing tsv file') }
    write.table(out, output.path, sep='\t', quote=F)
  }
}
