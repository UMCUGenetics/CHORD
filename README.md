# CHORD: Classifier of Homologous Recombination Deficiency

## Description
CHORD is a random forest model that uses the relative counts of somatic mutation contexts to predict
homologous recombination deficiency (HRD). The primary contexts used by CHORD are deletions with
flanking microhomology and 1-100kb structural duplications. Additionally, 1-100kb structural
duplications are used to distinguish BRCA1-type HRD from BRCA2-type HRD.

## Installation
Use devtools to install CHORD and mutSigExtractor directly from github. mutSigExtractor is required
for extracting the features that CHORD uses.
```
install.packages("devtools")
library(devtools)

install_github('https://github.com/luannnguyen/mutSigExtractor')
install_github('https://github.com/luannnguyen/CHORD/')
```

## Usage and tutorial
Ideally, the inputs for prediction are vcf files containing **somatic (no germline)** SNVs, indels,
and SVs per sample.

However, it is also possible to run CHORD on non-standard vcfs or from other sources. To do this we
can create dataframes that ```extractSigsSnv()```, ```extractSigsIndel()```, and
```extractSigsSv()``` accept (functions from the ```mutSigExtractor``` package). The
output from these functions can then be gathered and passed to ```chordPredict()```.

To get started with CHORD, please see the [tutorial](http://htmlpreview.github.io/?https://github.com/luannnguyen/CHORD/blob/master/example/run_chord.html).



