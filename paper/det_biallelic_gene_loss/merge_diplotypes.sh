#!/bin/bash

base_dir=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/
data_dir=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD_data/
vcf_subset_dir=$data_dir/HMF_DR010_DR047/vcf_subset/

out_txt=$base_dir/HMF_DR010_DR047/scripts/annotate_genes/hmf_gene_diplotypes_max.txt.gz

#touch $out_txt
counter=0
for i in $vcf_subset_dir/*; do
	counter=$((counter+1))

	in_txt=$i/gene_statuses/gene_diplotypes_max.txt.gz
	sample_name=$(basename $i)
	echo -ne "Processing [$counter]: $sample_name\r"

	## Write header using first file
	if [[ $counter -eq 1 ]]; then
		zcat $in_txt | head -n 1 | 
		awk '{print "sample""\t"$0}' | 
		gzip -c > $out_txt
	fi

	zcat $in_txt | tail -n +2 | 
	#grep BRCA | 
	awk -v sample_name="$sample_name" '{print sample_name"\t"$0}' | 
	gzip -c >> $out_txt

done

echo -e '\n'
