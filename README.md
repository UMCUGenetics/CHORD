# CHORD: Classifier of Homologous Recombination Deficiency

## Description
CHORD is a random forest model that uses the relative counts of somatic mutation contexts to predict
homologous recombination deficiency (HRD). The primary contexts used by CHORD are deletions with
flanking microhomology and 1-100kb structural duplications. Additionally, 1-100kb structural
duplications are used to distinguish BRCA1-type HRD from BRCA2-type HRD.

## Installation
CHORD is itself an R package. Before using CHORD, some other R packages will also need to be 
installed, with the main one being mutSigExtractor, which is required for extracting the mutation 
contexts that CHORD uses. The below code can be used to install CHORD and mutSigExtractor, as well 
as the dependencies for both packages.

```
## Bioconductor packages required by mutSigExtractor
install.packages('BiocManager')
BiocManager::install(
	c('GenomicRanges','VariantAnnotation','BSgenome','BSgenome.Hsapiens.UCSC.hg19')
)

## randomForest is required by CHORD
install.packages('randomForest')

## Install CHORD and mutSigExtractor directly from github using devtools
install.packages("devtools")
library(devtools)

install_github('https://github.com/luannnguyen/mutSigExtractor')
install_github('https://github.com/luannnguyen/CHORD/')
```

## Usage and tutorial
Ideally, the inputs for using CHORD are vcf files containing **somatic (no germline)** SNVs, indels,
and SVs per sample. ```extractSigsChord()``` is used to extract the relevant features (mutation
contexts) which are then passed to ```chordPredict()``` to make the HRD prediction.

However, it is also possible to run CHORD on non-standard vcfs or from other sources. To do this we
can create dataframes that ```extractSigsSnv()```, ```extractSigsIndel()```, and
```extractSigsSv()``` accept (functions from the ```mutSigExtractor``` package). The
output from these functions can then be gathered and passed to ```chordPredict()```.

To get started with CHORD, please see the [tutorial](http://htmlpreview.github.io/?https://github.com/luannnguyen/CHORD/blob/master/example/run_chord.html).



