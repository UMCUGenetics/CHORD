# CHORD: Classifier of Homologous Recombination Deficiency

CHORD is a random forest model that uses various mutational signatures to predict homologous
recombination deficiency (HRD). Per sample, the required inputs for prediction are two vcf files,
one containing SNVs and indels, and one containing SVs. Ideally, these vcfs will have been produced
using the HMF variant calling pipeline (https://github.com/hartwigmedical/pipeline/tree/v4.8), which uses
Strelka for somatic SNV/indel calling, and GRIDSS for SV calling.

The primary feature used by CHORD is deletions with flanking microhomology as well as structural 
duplications (1-10kb and 10-100kb in length). Additionally, structural duplications are used to 
distinguish BRCA1-type HRD from BRCA2-type HRD.

The CHORD package is dependent on the R packages mutSigExtractor (
https://github.com/luannnguyen/mutSigExtractor) and randomForest, so be sure to
have these installed. The below code can be used to install these dependencies locally.
```
install.packages("randomForest")

## Use devtools to install mutSigExtractor directly from github
install.packages("devtools"); library(devtools)
install_github('https://github.com/luannnguyen/mutSigExtractor')
```

Predicting HRD is performed in 2 steps. First, the counts of specific mutation contexts are extracted from 
vcf (or compressed vcf.gz) files. Note that with many samples/large vcfs, it might be a good idea to run 
this step on an HPC.
```
## extractSigsChord() will extract the mutation contexts for one sample. This will return a one row dataframe.
contexts <- extractSigsChord(
   vcf.snv = '/path/to/vcf_with_snvs',
   vcf.indel = '/path/to/vcf_with_indels',
   vcf.sv = '/path/to/vcf_with_svs',
   sample.name = 'can_be_anything',
   output.path = NULL ## If this is not NULL, the output dataframe will be written to this path
)
```

Then, the mutation context counts are run through the model to make the prediction.
```
## With multiple samples, the one row dataframes outputted by extractSigsChord() can be merged into 
## one dataframe. For example:
l_contexts <- list(contexts1, contexts2, contexts3) ## In reality this list will probably be generated from an lapply loop
contexts <- do.call(rbind, l_contexts)

## Prediction
pred <- chordPredict(contexts)
```