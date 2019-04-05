# CHORD: Classifier of Homologous Recombination Deficiency

CHORD is a random forest model that uses various mutational signatures to predict homologous
recombination deficiency (HRD). Per sample, the required inputs for prediction are two vcf files,
one containing SNVs and indels, and one containing SVs. Ideally, these vcfs will have been produced
using the HMF variant calling pipeline v4.4 (https://github.com/hartwigmedical/pipeline), which uses
Strelka for somatic SNV/indel calling, and Manta + BPI (Breakpoint Inspector; custom code for SV
post-processing).

The primary feature used by CHORD is deletions with flanking microhomology. Also used structural 
duplications (1,000-10,000bp & 10,000-100,000bp in length). Structural duplications are used to 
distinguish BRCA1(-like) from BRCA2(-like) HRD. More info about CHORD (training, performance, etc...) 
can be found in info/chord_training_and_perf.pptx.

The CHORD package is dependent on the R packages mutSigExtractor and randomForest, so be sure to
have these installed. The below code can be used to install these dependencies locally.
```
install.packages("randomForest")

## Use devtools to install mutSigExtractor directly from github
install.packages("devtools"); library(devtools)
install_github('https://github.com/luannnguyen/mutSigExtractor')
```

Predicting HRD is performed in 2 steps. First, signatures are extracted from vcf (or compressed 
vcf.gz) files. Note that with many samples/large vcfs, it might be a good idea to run this step on 
an HPC.
```
## extractSigsChord() will extract signatures for one sample. This will return a one row dataframe.
sigs <- extractSigsChord(
   vcf.snv = '/path/to/vcf_with_snvs',
   vcf.indel = '/path/to/vcf_with_indels',
   vcf.sv = '/path/to/vcf_with_svs',
   sample.name = 'can_be_anything',
   output.path = NULL ## If this is not NULL, the output dataframe will be written to this path
)
```

Then, the signature values are run through the model to make the prediction.
```
## With multiple samples, the one row dataframes outputted by extractSigsChord() can be merged into 
## one dataframe. For example:
df_list <- list(df1, df2, df3) ## In reality this list will probably be generated from an lapply loop
sigs <- do.call(rbind, df_list)

## Prediction
pred <- chordPredict(sigs)
```