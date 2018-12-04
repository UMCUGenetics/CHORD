# CHORD (Classifier for Homologous Recombination Deficiency)

CHORD is a random forest model that was trained on mutational signatures to predict homologous
recombination deficiency (HRD). The main features used by the model are deletions with flanking
mirochomology, and the SV signatures 1,3 and 5 (described here: 
https://media.nature.com/original/nature-assets/nature/journal/v534/n7605/extref/nature17676-s3.zip, 
in Supplementary.Table.21.Signatures.v3.xlsx)

The CHORD package is dependent on the R packages mutSigExtractor and randomForest, so be sure to
have these installed.

Predicting HRD is performed in 2 steps.

First, signatures are extracted from vcf (or compressed vcf.gz) files.
```
sigs <- extractSigsChord(
   vcf.snv = '/path/to/vcf_with_snvs',
   vcf.indel = '/path/to/vcf_with_indels',
   vcf.sv = '/path/to/vcf_with_svs',
   sample.name = 'can_be_anything'
)
```

Then, the signature values are run through the model to make the prediction.
```
pred <- chordPredict(sigs)
```