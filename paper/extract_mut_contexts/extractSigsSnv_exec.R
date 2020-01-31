#!/usr/bin/env Rscript
## Run on hpc
library(devtools)
load_all('/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/mutSigExtractor/')

## HPC
args <- commandArgs(trailingOnly = TRUE)

#message('Extracting snv signatures...')
sigs_snv <- extractSigsSnv(vcf.file = args[1],
                           sample.name = args[2],
                           vcf.filter = 'PASS',
                           output = 'contexts')

#message('Writing to file...')
write.table(sigs_snv, args[3], sep = '\t', quote = F)



