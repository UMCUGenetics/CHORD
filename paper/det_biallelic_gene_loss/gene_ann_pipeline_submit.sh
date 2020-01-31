#!/bin/bash

base_dir=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/
data_dir=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD_data/

pipeline_sh=$base_dir/scripts_main/hmfGeneAnnotation/scripts/pipeline/pipeline.sh; source $pipeline_sh
bed_path=$base_dir/scripts_main/hmfGeneAnnotation/data/gene_selection/genes.bed
dgs_ini_path=$base_dir/scripts_main/hmfGeneAnnotation/scripts/pipeline/detGeneStatuses_ini.R

#========= Submit manifest =========#
hmf_data_dir=/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-DR047/data/

manifest_path=$base_dir/HMF_DR010_DR047/manifest/manifest.txt
variants_dir=$data_dir/HMF_DR010_DR047/vcf_subset/; mkdir -p $variants_dir

counter=0
tail $manifest_path -n +2 | while read sample_name set_name germ_vcf som_vcf sv_vcf gene_cnv purity; do
	counter=$((counter+1))

	#if [[ $counter -le 2700 ]]; then continue; fi
	
	echo -e "\n### [$counter] Submitting pipeline for $sample_name ###"
	
	out_dir=$variants_dir/$sample_name; mkdir -p $out_dir

	#--------- inputs ---------#
	bed_path=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/data/gene_selection/genes.bed
	#purity_path=$hmf_data_dir/$set_name/$purity
	gene_cnv_path=$hmf_data_dir/$set_name/$gene_cnv

	germ_vcf_path=$hmf_data_dir/$set_name/$germ_vcf
	som_vcf_path=$hmf_data_dir/$set_name/$som_vcf

	#--------- submit ---------#
	pipeline -o $out_dir -b $bed_path -n $sample_name \
	-c $gene_cnv_path -g $germ_vcf_path -s $som_vcf_path \
	-i $dgs_ini_path -t 1 -k 7

	#if [[ $counter -eq 1 ]]; then break; fi
done

