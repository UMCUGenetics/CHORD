CHORD: Classifier of HOmologous Recombination Deficiency
================

## Description

CHORD is a random forest model that uses the relative counts of somatic
mutation contexts to predict homologous recombination deficiency (HRD).
The primary contexts used by CHORD are deletions with flanking
microhomology and 1-100kb structural duplications. Additionally, 1-100kb
structural duplications are used to distinguish BRCA1-type HRD from
BRCA2-type HRD.

## Installation

CHORD is itself an R package. Before using CHORD, some other R packages
will also need to be installed, with the main one being mutSigExtractor,
which is required for extracting the mutation contexts that CHORD uses.
The below code can be used to install CHORD and mutSigExtractor, as well
as the dependencies for both packages.

``` r
## Bioconductor packages required by mutSigExtractor
install.packages('BiocManager')
BiocManager::install('BSgenome') ## Install genome parser
BiocManager::install('BSgenome.Hsapiens.UCSC.hg19') ## Install relevant genome

## randomForest is required by CHORD
install.packages('randomForest')

## Install CHORD and mutSigExtractor directly from github using devtools
install.packages("devtools")
library(devtools)

install_github('https://github.com/luannnguyen/mutSigExtractor')
install_github('https://github.com/luannnguyen/CHORD/')
```

## Quick start

### Vcfs as input

CHORD can be supplied with vcf files containing **somatic (no
germline)** SNVs, indels, and SVs per sample. The first step is to
extract the mutation contexts required by CHORD:

``` r
contexts <- extractSigsChord(
  vcf.snv = '/path/to/vcf/with/snvs/',
  vcf.indel = '/path/to/vcf/with/snvs/',
  vcf.sv = '/path/to/vcf/with/svs/',
  sv.caller = 'gridss' ## Can be 'gridss' (default) or 'manta'
)
```

Note that there is no real standard for SV vcfs. However, parsing vcfs
produced by Manta or GRIDSS has been built into CHORD.

### Dataframes as input

Alternatively, dataframes can be used as input (handy if you want to
parse the vcfs yourself and supply this preprocessed data to CHORD).

  - For SNVs and indels, a dataframe with the columns: chrom, pos, ref,
    alt
  - For SVs, a dataframe with the columns: sv\_type, sv\_len. For
    sv\_type, values must be DEL, DUP, INV, TRA

Note that the column names themselves do not matter, as long as the
columns are in the aforementioned order.

``` r
contexts <- extractSigsChord(
  df.snv = bed_file_like_dataframe_with_snvs,
  df.indel = bed_file_like_dataframe_with_indels,
  df.sv = dataframe_with_svtype_and_sv_len
)
```

### Inputs are flexible

Often, SNV and indels are reported in the same vcfs. In such cases, the
vcf path can be specified to the `vcf.snv` argument (`vcf.indel` can be
left out).

``` r
contexts <- extractSigsChord(
  vcf.snv = '/path/to/vcf/with/snvs_and_indels/',
  vcf.sv = '/path/to/vcf/with/svs/',
  sv.caller = 'gridss' ## Can be 'gridss' (default) or 'manta'
)
```

The vcf.\* and df.\* arguments can be mixed. This is especially useful
in the case of SV vcfs, since it is possible that you will need to parse
these yourself.

``` r
contexts <- extractSigsChord(
  vcf.snv = '/path/to/vcf/with/snvs/',
  vcf.indel = '/path/to/vcf/with/snvs/',
  df.sv = dataframe_with_svtype_and_sv_len
)
```

A different reference genome than the default can be used. Genomes
should be BSgenomes.

``` r
## Make sure to install and load the desired first
install.packages('BiocManager')
ref_genome <- 'BSgenome.Hsapiens.UCSC.hg38'
BiocManager::install(ref_genome)
library(ref_genome)

## Specify the name of the BSgenome object to the ref.genome argument
contexts <- extractSigsChord(
  vcf.snv = '/path/to/vcf/with/snvs/',
  vcf.indel = '/path/to/vcf/with/snvs/',
  vcf.sv = '/path/to/vcf/with/svs/',
  sv.caller = 'gridss' ## Can be 'gridss' (default) or 'manta',
  ref.genome=ref_genome
)
```

### Predicting HRD

Once we have the mutation contexts, a prediction can be made.

    chordPredict(contexts)

For more details on how to use CHORD, please see the below tutorial.

## Tutorial

This tutorial will demonstrate how HRD prediction can be performed
locally. Somatic vcfs from a few primary tumors from the 560 breast
cancer cohort (BRCA-EU, ICGC) will be used as example input. These are
located at `example/vcf/` in the CHORD git repository. For large vcfs
and/or many samples however, it is advised to run CHORD on a
high-performance cluster (HPC).

In part 1, we will use the vcfs directly as input for CHORD.

Part 2 of this tutorial will provide a demonstration of how we can parse
the vcfs into dataframes that CHORD accepts. This is particularly useful
if you may have non-standard vcfs, as is often the case with SV vcfs. Of
course, it is possible to produce the required input dataframes from any
data source, and not just vcfs; it only requires a bit of data wrangling
:).

To run the example, set the working directory to the `example/`
directory in the CHORD package.

``` r
setwd('./CHORD/example/')
```

### 1\. Running CHORD on vcfs

#### Loading vcfs into R

We will first create a dataframe to store the vcf paths as well as
assign them sample names.

``` r
library(CHORD)

options(stringsAsFactors=F)

## Get file paths from the vcf/ subdirectory
vcf_files <- data.frame(
  snv_indel=list.files('vcf', pattern='*snv_indel.vcf.gz', full.names=T),
  sv=list.files('vcf', pattern='*sv.vcf.gz', full.names=T)
)

## Assign sample names to files
vcf_files$sample <- sapply(strsplit(basename(vcf_files$snv_indel),'_'),`[`,1)

vcf_files
```

    ##                      snv_indel                    sv  sample
    ## 1 vcf/PD10010_snv_indel.vcf.gz vcf/PD10010_sv.vcf.gz PD10010
    ## 2 vcf/PD11352_snv_indel.vcf.gz vcf/PD11352_sv.vcf.gz PD11352
    ## 3  vcf/PD3905_snv_indel.vcf.gz  vcf/PD3905_sv.vcf.gz  PD3905
    ## 4  vcf/PD4116_snv_indel.vcf.gz  vcf/PD4116_sv.vcf.gz  PD4116
    ## 5  vcf/PD7344_snv_indel.vcf.gz  vcf/PD7344_sv.vcf.gz  PD7344

Vcfs can either be in plain text format (.vcf) or compressed (.gz).

SNVs and indels are often stored in the same vcf (as is the case with
the example vcfs). CHORD therefore accepts such vcfs as input; SNVs and
indels are automatically split up. However, it is also possible to
specify the SNV and indel vcfs separately (by specifying the paths
separately to `vcf.snv` and `vcf.indel` to `extractSigsChord()`; see
below).

#### Extracting mutation contexts

Create a directory to write all subsequent output.

``` r
dir.create('output/')
```

***Note***: the output from this tutorial exists in the
`example/output/` directory by default. Delete this directory to create
new output files.

Extract the features that are used as input for CHORD. With the example
data, it should take about 30s to run.

``` r
## Make dir to output contexts
contexts_dir <- 'output/contexts/'
dir.create(contexts_dir, recursive=T)

## Extract contexts for all samples
for(i in 1:nrow(vcf_files)){
  
  params <- as.list(vcf_files[i,])
  out_path <- paste0(contexts_dir,'/',params$sample,'_contexts.txt')
  
  if(!file.exists(out_path)){
    extractSigsChord(
      vcf.snv=params$snv_indel,
      vcf.sv=params$sv, sv.caller='manta',
      sample.name=params$sample,
      output.path=out_path, verbose=F
    )
  }
}
```

The `sv.caller` parameter has been included in `extractSigsChord()`
since different SV callers report SVs in different ways. In CHORD, we
have included parsing for vcf produced by Manta and GRIDSS
(`sv.caller='manta'` or `sv.caller='gridss'`).

We can then merge the contexts into a matrix.

``` r
## Get the paths of the contexts txt files
context_files <- list.files(contexts_dir, full.names=T)

## Read the txt files into R
l_contexts <- lapply(context_files, function(i){
  read.delim(i, check.names=F)
})

## Merged the contexts into a matrix
merged_contexts <- do.call(rbind, l_contexts)

## Write to output directory
write.table(merged_contexts, 'output/merged_contexts.txt', sep='\t', quote=F)
```

Below we can see what a part of the context matrix looks like.

``` r
merged_contexts[,1:5]
```

    ##         A[C>A]A A[C>A]C A[C>A]G A[C>A]T C[C>A]A
    ## PD10010      37      25       8      24      34
    ## PD11352      29      40       3      21      30
    ## PD3905       89      61      14      61      91
    ## PD4116      123     136      15     142     164
    ## PD7344       10       9       1       6      11

#### Predicting HRD and interpreting CHORD’s output

Once we have the context matrix ready, we can use it for predicting HRD.

``` r
pred <- chordPredict(merged_contexts, verbose=F)
write.table(pred, 'output/chord_pred.txt', sep='\t', quote=F)

pred
```

    ##    sample p_none p_BRCA1 p_BRCA2 p_hrd             qc is_hrd   hrd_type
    ## 1 PD10010  1.000   0.000   0.000 0.000           pass  FALSE       none
    ## 2 PD11352  1.000   0.000   0.000 0.000    low_sv_load  FALSE       <NA>
    ## 3  PD3905  0.060   0.924   0.016 0.940           pass   TRUE BRCA1_type
    ## 4  PD4116  0.102   0.106   0.792 0.898           pass   TRUE BRCA2_type
    ## 5  PD7344  0.998   0.000   0.002 0.002 low_indel_load     NA       <NA>

#### Interpreting the output

##### Main output

CHORD outputs the probability of:

  - HRD: `p_hrd`
  - BRCA1-type HRD: `p_BRCA1`
  - BRCA2-type HRD: `p_BRCA2`
  - Not HRD: `p_none`

Note that `p_hrd` = `p_BRCA1` + `p_BRCA2`. Furthermore, `p_none` +
`p_BRCA1` + `p_BRCA2` = 1.

`is_hrd` tells us if `p_hrd` **\>= 0.5 (the HRD classification
threshold)**. If this is TRUE, then `hrd_type` will tell us if the
sample has BRCA1-type HRD or BRCA2-type HRD (=max(`p_BRCA1`,`p_BRCA2`)).

##### QC checks

Under `qc`, we can see that PD7344 has `low_indel_load`. CHORD requires
\>=50 indels to accurately determine whether a sample is HRD. PD7344
fails this check and therefore `is_hrd` is NA. Consequently, HRD subtype
can not be determined, and thus `hrd_type` is NA.

Also under `qc`, we can see that PD11352 has `low_sv_load`. CHORD
requires \>=30 SVs to accurately determine HRD subtype. Thus, while this
sample had sufficient number of indels to determine HRD, the number of
SVs was not sufficient, and therefore `hrd_type` is NA.

The user may of course ignore these QC warnings and proceed with using
the raw probabilities outputted by CHORD (`p_hrd` and/or
`p_BRCA1`/`pBRCA2`).

### 2\. Running CHORD from dataframes

To run CHORD on non-standard vcfs or from other sources, we can create
dataframes that `extractSigsChord()` accepts. The output can then be
passed to `chordPredict()`.

In this part of the tutorial, we will create the required dataframes
from the example vcfs. This will thus also serve as a short
demonstration on how one could parse vcfs in R. Furthermore, we will
focus on one sample for simplicity’s sake.

#### SNV and indel dataframes

For SNVs and indels, a dataframe containing CHROM, POS, REF and ALT
information (i.e. a bed file like format) is required. The column names
of this dataframe can be anything; as long as the columns are in the
aforementioned order.

We can use `readVcfFields()` from the `mutSigExtractor` package to read
a vcf into R as a
dataframe.

``` r
df_snv_indel <- readVcfFields('vcf/PD3905_snv_indel.vcf.gz', fields=c('CHROM','POS','REF','ALT'))
```

Alternatively, we can specify the vcf columns we want by
index.

``` r
df_snv_indel <- readVcfFields('vcf/PD3905_snv_indel.vcf.gz', fields=c(1,2,4,5))
```

The dataframe shown below is in the correct format as required by CHORD.
This dataframe contains both SNVs and indels. This is not a problem as
we can specify the argument`df.snv` (and leaving out `df.indel`) in
`extractSigsChord()` as is done later in the tutorial.

``` r
df_snv_indel[25:30,]
```

    ##    CHROM      POS           REF ALT
    ## 25     1 16145827             A   C
    ## 26     1 16333830 GTGTGAAGGTGTC   G
    ## 27     1 16492085             G   C
    ## 28     1 17890303             C   G
    ## 29     1 18877885             G   A
    ## 30     1 18919776             T   C

#### Extracting SV contexts

For SVs, a dataframe is required with SV type (values must be one of the
following: DEL, DUP, INV, TRA) and SV length information. Again, the
column names of this dataframe can be anything; as long as the columns
are in the aforementioned order.

Again, `readVcfFields()` can be used to read the vcf. The example vcfs
have been produced by Manta, which reports SVTYPE and SVLEN in the INFO
field, so only this field needs to be extracted from the vcf.

``` r
vcf_sv <- readVcfFields('vcf/PD3905_sv.vcf.gz', fields='INFO')
print(vcf_sv[1:5,,drop=F], right=F)
```

    ##   INFO                                                                                                                                                 
    ## 1 SVLEN=-13737754;SVTYPE=TRA;CHR2=3;END=10206978;Donor_ID=DO224557;project=BRCA-EU;icgc_sample_ID=PD3905a;icgc_specimen_ID=SA569900;assembly_version=NA
    ## 2 SVLEN=1696;SVTYPE=DEL;CHR2=1;END=54899207;Donor_ID=DO224557;project=BRCA-EU;icgc_sample_ID=PD3905a;icgc_specimen_ID=SA569900;assembly_version=NA     
    ## 3 SVLEN=22644;SVTYPE=DEL;CHR2=1;END=72263583;Donor_ID=DO224557;project=BRCA-EU;icgc_sample_ID=PD3905a;icgc_specimen_ID=SA569900;assembly_version=NA    
    ## 4 SVLEN=1703;SVTYPE=DUP;CHR2=1;END=75618395;Donor_ID=DO224557;project=BRCA-EU;icgc_sample_ID=PD3905a;icgc_specimen_ID=SA569900;assembly_version=NA     
    ## 5 SVLEN=1789;SVTYPE=DEL;CHR2=1;END=75619297;Donor_ID=DO224557;project=BRCA-EU;icgc_sample_ID=PD3905a;icgc_specimen_ID=SA569900;assembly_version=NA

Then, split the data each INFO entry, which are separated by `;`.

``` r
vcf_sv_info <- strsplit(vcf_sv$INFO,';')
vcf_sv_info[1:2]
```

    ## [[1]]
    ## [1] "SVLEN=-13737754"           "SVTYPE=TRA"               
    ## [3] "CHR2=3"                    "END=10206978"             
    ## [5] "Donor_ID=DO224557"         "project=BRCA-EU"          
    ## [7] "icgc_sample_ID=PD3905a"    "icgc_specimen_ID=SA569900"
    ## [9] "assembly_version=NA"      
    ## 
    ## [[2]]
    ## [1] "SVLEN=1696"                "SVTYPE=DEL"               
    ## [3] "CHR2=1"                    "END=54899207"             
    ## [5] "Donor_ID=DO224557"         "project=BRCA-EU"          
    ## [7] "icgc_sample_ID=PD3905a"    "icgc_specimen_ID=SA569900"
    ## [9] "assembly_version=NA"

Extract the SVTYPE and SVLEN data, and put it in the correct dataframe
format.

``` r
sv_type <- sapply(vcf_sv_info,`[`,2) ## Select the 2nd object from each INFO entry
sv_type <- gsub('SVTYPE=','',sv_type) ## Remove the 'SVTYPE=' prefix

sv_len <- sapply(vcf_sv_info,`[`,1) ## Select the 1st object from each INFO entry
sv_len <- gsub('SVLEN=','',sv_len) ## Remove the 'SVLEN=' prefix
sv_len <- as.integer(sv_len) ## Convert the character vector to an integer vector

df_sv <- data.frame(sv_type, sv_len)
head(df_sv)
```

    ##   sv_type    sv_len
    ## 1     TRA -13737754
    ## 2     DEL      1696
    ## 3     DEL     22644
    ## 4     DUP      1703
    ## 5     DEL      1789
    ## 6     DEL     49256

You may notice that translocations (TRA) have an SV length, which
doesn’t really make sense. This is no problem since SV length info is
discarded for translocations.

#### Predicting HRD

Extract the relevant mutation contexts from the variant data.

``` r
contexts <- extractSigsChord(
  df.snv = df_snv_indel,
  df.sv = df_sv,
  sample.name='PD3905'
)
```

Below we can see what a part of the context matrix looks
    like.

``` r
contexts[,c(1:4, 97:104, 127:131),drop=F]
```

    ##        A[C>A]A A[C>A]C A[C>A]G A[C>A]T del.rep.len.1 del.rep.len.2
    ## PD3905      89      61      14      61            24             6
    ##        del.rep.len.3 del.rep.len.4 del.rep.len.5 ins.rep.len.1 ins.rep.len.2
    ## PD3905             4             1             0            19             1
    ##        ins.rep.len.3 DEL_0e00_1e03_bp DEL_1e03_1e04_bp DEL_1e04_1e05_bp
    ## PD3905             0                0                8                8
    ##        DEL_1e05_1e06_bp DEL_1e06_1e07_bp
    ## PD3905                1                2

Lastly, make the HRD prediction.

``` r
pred <- chordPredict(contexts)
pred
```

    ##   sample p_none p_BRCA1 p_BRCA2 p_hrd   qc is_hrd   hrd_type
    ## 1 PD3905   0.06   0.924   0.016  0.94 pass   TRUE BRCA1_type

Please refer back to section 1 of the tutorial for interpreting CHORD’s
output.
