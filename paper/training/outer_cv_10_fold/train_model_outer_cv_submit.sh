#!/bin/bash

## qlogin -l tmpspace=50G -l h_rt=10:00:00 -pe threaded 20

base_dir=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/

script_path=$base_dir/HMF_DR010_DR047/scripts/chord_training/main/scripts/train_model_outer_cv.R

submitTraining (){
	fold_dir=$1
	guixr load-profile ~/.guix-profile << EOF1
Rscript $script_path $fold_dir
EOF1
}

parent_dir=$base_dir/HMF_DR010_DR047/scripts/chord_training/main/output/60_CV_snvContext_svContext_allSigsRel_noSuppEvidence_noBoruta/
for i in $parent_dir/fold_*; do
	submitTraining $i
done

