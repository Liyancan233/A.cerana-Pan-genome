#!/usr/local/bin/bash
#conda base
#in sample_results,set up dif floder by reads len
#sample 4X 5X 6X apiscerana reads
sample_result_dir=${1}

cd ${sample_result_dir}/200BP

reads_list_one=$(find ${sample_result_dir}/200BP/*_1.fastq.gz | sed 's#.*/##' | awk 'BEGIN{FS="_";OFS="_"}{print $1,$2,$3}')
for x in ${reads_list_one};
do
    nohup seqkit sample ${x}_1.fastq.gz -n 3300000 -s 11 -o 3x_${x}_1.fq.gz &
    nohup seqkit sample ${x}_2.fastq.gz -n 3300000 -s 11 -o 3x_${x}_2.fq.gz &
    nohup seqkit sample ${x}_1.fastq.gz -n 4400000 -s 11 -o 4x_${x}_1.fq.gz &
    nohup seqkit sample ${x}_2.fastq.gz -n 4400000 -s 11 -o 4x_${x}_2.fq.gz &
    nohup seqkit sample ${x}_1.fastq.gz -n 5100000 -s 11 -o 5x_${x}_1.fq.gz &
    nohup seqkit sample ${x}_2.fastq.gz -n 5100000 -s 11 -o 5x_${x}_2.fq.gz &
    nohup seqkit sample ${x}_1.fastq.gz -n 6200000 -s 11 -o 6x_${x}_1.fq.gz &
    nohup seqkit sample ${x}_2.fastq.gz -n 6200000 -s 11 -o 6x_${x}_2.fq.gz &
done


echo " 200BP is done !!!"

cd ${sample_result_dir}/250BP

reads_list_two=$(find ${sample_result_dir}/250BP/*_1.fastq.gz | sed 's#.*/##' | awk 'BEGIN{FS="_";OFS="_"}{print $1,$2,$3}')
for x in ${reads_list_two};
do
    nohup seqkit sample ${x}_1.fastq.gz -n 2640000 -s 11 -o 3x_${x}_1.fq.gz &
    nohup seqkit sample ${x}_2.fastq.gz -n 2640000 -s 11 -o 3x_${x}_2.fq.gz &
    nohup seqkit sample ${x}_1.fastq.gz -n 3520000 -s 11 -o 4x_${x}_1.fq.gz &
    nohup seqkit sample ${x}_2.fastq.gz -n 3520000 -s 11 -o 4x_${x}_2.fq.gz &
    nohup seqkit sample ${x}_1.fastq.gz -n 4080000 -s 11 -o 5x_${x}_1.fq.gz &
    nohup seqkit sample ${x}_2.fastq.gz -n 4080000 -s 11 -o 5x_${x}_2.fq.gz &
    nohup seqkit sample ${x}_1.fastq.gz -n 4960000 -s 11 -o 6x_${x}_1.fq.gz &
    nohup seqkit sample ${x}_2.fastq.gz -n 4960000 -s 11 -o 6x_${x}_2.fq.gz &
done


echo " 250BP is done !!!"

cd ${sample_result_dir}/300BP

reads_list_three=$(find ${sample_result_dir}/300BP/*_1.fastq.gz | sed 's#.*/##' | awk 'BEGIN{FS="_";OFS="_"}{print $1,$2,$3}')
for x in ${reads_list_three};
do
    nohup seqkit sample ${x}_1.fastq.gz -n 2200000 -s 11 -o 3x_${x}_1.fq.gz &
    nohup seqkit sample ${x}_2.fastq.gz -n 2200000 -s 11 -o 3x_${x}_2.fq.gz &
    nohup seqkit sample ${x}_1.fastq.gz -n 2930000 -s 11 -o 4x_${x}_1.fq.gz &
    nohup seqkit sample ${x}_2.fastq.gz -n 2930000 -s 11 -o 4x_${x}_2.fq.gz &
    nohup seqkit sample ${x}_1.fastq.gz -n 3400000 -s 11 -o 5x_${x}_1.fq.gz &
    nohup seqkit sample ${x}_2.fastq.gz -n 3400000 -s 11 -o 5x_${x}_2.fq.gz &
    nohup seqkit sample ${x}_1.fastq.gz -n 4130000 -s 11 -o 6x_${x}_1.fq.gz &
    nohup seqkit sample ${x}_2.fastq.gz -n 4130000 -s 11 -o 6x_${x}_2.fq.gz &
done

















