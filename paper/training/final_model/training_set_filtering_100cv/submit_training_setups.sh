#!/bin/bash

## qlogin -l tmpspace=50G -l h_rt=10:00:00 -pe threaded 20

base_dir=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/
parent_dir=$base_dir//HMF_DR010_DR047/scripts/chord_training/main/output/60_snvContext_svContext_allSigsRel_noSuppEvidence_customBlacklist3_noBoruta/cv100/replicates/

job_dir=$parent_dir/jobs/; mkdir -p $job_dir; cd $job_dir

for i in $parent_dir/replicates/*/train_chord_custom.R; do

	job_file=$job_dir/$(basename $(dirname $i)).job
	done_file=$(dirname $i)/done

echo "guixr load-profile ~/.guix-profile --<<EOF
Rscript $i && touch $done_file
EOF" > $job_file
	
	if [[ ! -f $done_file ]]; then
		qsub -cwd -l h_rt=0:30:00 -l h_vmem=8G -pe threaded 12 $job_file
	fi
done
