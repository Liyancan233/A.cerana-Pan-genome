#!/usr/local/bin/bash
ref_dir=${1}  #e.g./mnt/cheng/lyc/genome/TMP
out_put_dir=${2} #e.g./mnt/cheng/lyc/genome
#genomelist=$(ls ${ref_dir}/* | awk 'BEGIN{FS="_";OFS="_"}{print $0}')

#path=$ref_dir
files=$(ls $ref_dir)
for x in $files;
do
    cd $ref_dir/${x}
    mv ${x}_output.g.vcf.gz $out_put_dir
    mv ${x}_output.g.vcf.gz.tbi $out_put_dir
done

