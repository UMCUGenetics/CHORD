#!/usr/bin/env Rscript
## Run on hpc
library(devtools)
load_all('/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/mutSigExtractor/')

## HPC
extractSigsSv_exec <- function(vcf.sv, vcfname.sv, out.path.sv){

   sigs_sv <- extractSigsSv(vcf.sv, vcf.filter='PASS', sv.caller='gridss',
                            sample.name = vcfname.sv, verbose = F, output = 'contexts',
                            sv.len.cutoffs = c(0, 10^3, 10^4, 10^5, 10^6, 10^7,Inf))
   write.table(sigs_sv, out.path.sv, sep = '\t', quote = F)

}

args <- commandArgs(trailingOnly = TRUE)

extractSigsSv_exec(args[1], args[2], args[3])