# CHORD (Classifier for Homologous Recombination Deficiency)

CHORD is a random forest model that was trained on mutational signatures to predict homologous
recombination deficiency (HRD). The main features used by the model are deletions with flanking
mirochomology, and the SV signatures 1,3 and 5 (described here: 
https://media.nature.com/original/nature-assets/nature/journal/v534/n7605/extref/nature17676-s3.zip, 
in Supplementary.Table.21.Signatures.v3.xlsx)

The CHORD package is dependent on the R packages mutSigExtractor and randomForest, so be sure to
have these installed.

Predicting HRD is performed in 2 steps. First, signatures are extracted from vcf (or compressed vcf.gz) 
files. Note that with many samples/large vcfs, it might be a good idea to run this step on an HPC.
```
## extractSigsChord() will extract signatures for one sample. This will return a one row dataframe.
sigs <- extractSigsChord(
   vcf.snv = '/path/to/vcf_with_snvs',
   vcf.indel = '/path/to/vcf_with_indels',
   vcf.sv = '/path/to/vcf_with_svs',
   sample.name = 'can_be_anything'
)
```

Then, the signature values are run through the model to make the prediction.
```
## With multiple samples, the one row dataframes outputted by extractSigsChord() can be merged into one 
## dataframe. For example:
df_list <- list(output_df1, output_df2, output_df3) ## in reality this list will probably be generated from an lapply loop
sigs <- do.call(rbind, df_list)

## Prediction
pred <- chordPredict(sigs)
```