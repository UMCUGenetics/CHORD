#!/usr/bin/env Rscript
## Run on hpc
library(devtools)
load_all('/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/mutSigExtractor/')

## HPC
args <- commandArgs(trailingOnly = TRUE)

#message('Extracting indel signatures...')
sigs_indel <- extractSigsIndel(vcf.file = args[1],
                               sample.name = args[2],
                               vcf.filter = 'PASS')

#message('Writing to file...')
write.table(sigs_indel, args[3], sep = '\t', quote = F)



