#!/bin/bash

base_dir=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/

#========= Submit manifest =========#
hmf_data_dir=/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-DR047/data/
manifest_path=$base_dir/HMF_DR010_DR047/manifest/manifest.txt

scripts_dir=$base_dir/HMF_DR010_DR047/scripts/extract_sigs/
r_snv=$scripts_dir/extractSigsSnv_exec.R
r_indel=$scripts_dir/extractSigsIndel_exec.R
r_sv=$scripts_dir/extractSigsSv_exec.R

job_dir=$base_dir/HMF_DR010_DR047/jobs/extractSigs/
mkdir -p $job_dir
cd $job_dir

out_dir=$base_dir/HMF_DR010_DR047/matrices/
out_snv=$out_dir/snv/; mkdir -p $out_snv
out_indel=$out_dir/indel/; mkdir -p $out_indel
out_sv=$out_dir/sv/; mkdir -p $out_sv

log_dir=$out_dir/log/
mkdir -p $log_dir
log_snv=$log_dir/snv.log; touch $log_snv
log_indel=$log_dir/indel.log; touch $log_indel
log_sv=$log_dir/sv.log; touch $log_sv

counter=0
tail $manifest_path -n +2 | while read sample_name set_name germ_vcf_name som_vcf_name sv_vcf_name gene_cnv_name purity_name; do
	counter=$((counter+1))

	echo -e "\n########## [$counter] Submitting extract sigs for: $sample_name ##########"

	## Get relative paths for vcfs
	vcf_som=$hmf_data_dir/$set_name/$som_vcf_name
	vcf_sv=$hmf_data_dir/$set_name/$sv_vcf_name

	## Make jobs
	job_snv=$job_dir/snv_${sample_name}.job
	job_indel=$job_dir/indel_${sample_name}.job
	job_sv=$job_dir/sv_${sample_name}.job

	## Submit
	qsubDefault (){
		#qsub -S /bin/bash -cwd -m ea -l h_rt=1:00:00 -l h_vmem=2G $1 
		qsub -S /bin/bash -cwd -m ea -l h_rt=12:00:00 -l h_vmem=6G $1 
	}

	if [[ ! -f $out_snv/${sample_name}_snv.txt && -f $vcf_som ]]; then 
		echo guixr load-profile ~/.guix-profile '--<<EOF' > $job_snv
		echo Rscript $r_snv $vcf_som $sample_name $out_snv/${sample_name}_snv.txt >> $job_snv
		echo 'EOF' >> $job_snv
		qsubDefault $job_snv
	fi
	
	if [[ ! -f $out_indel/${sample_name}_indel.txt && -f $vcf_som ]]; then 
		echo guixr load-profile ~/.guix-profile '--<<EOF' > $job_indel
		echo Rscript $r_indel $vcf_som $sample_name $out_indel/${sample_name}_indel.txt >> $job_indel
		echo 'EOF' >> $job_indel
		qsubDefault $job_indel
	fi
	
	if [[ ! -f $out_sv/${sample_name}_sv.txt && -f $vcf_sv ]]; then
		echo guixr load-profile ~/.guix-profile '--<<EOF' > $job_sv
		echo Rscript $r_sv $vcf_sv $sample_name $out_sv/${sample_name}_sv.txt >> $job_sv
		echo 'EOF' >> $job_sv
		qsubDefault $job_sv
	fi

	#if [[ $counter -eq 1 ]]; then break; fi

done
