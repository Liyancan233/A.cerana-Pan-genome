#!/usr/local/bin/bash
genome_dir=${1} 
genomelist=$(find ${genome_dir}/*vcf.gz | sed 's#.*/##' | awk 'BEGIN{FS="_";OFS="_"}{print $0}')
for x in $genomelist;
do
	gatk IndexFeatureFile -I $x
done

