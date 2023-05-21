#!/usr/local/bin/bash
#conda base
raw_reads_dir=${1}
#step1 select sample from all reads with *1/2.fastq.gz end
cd ${raw_reads_dir}
mkdir sample_result
genomelist=$(ls *1.fastq.gz | sort --random-sort | awk '{if(NR%25 ==1) print $0}' | awk 'BEGIN{FS="_";OFS="_"}{print $1,$2,$3}')
for x in $genomelist;
do
    mv ${x}* sample_result/
done
#set up two floder 200BP and 250 BP
#
