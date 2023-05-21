#!/usr/local/bin/bash
ref_dir=${1}  
mkdir -p ${ref_dir}/repeat_result
mkdir -p ${ref_dir}/repeat_result/5_rep
mkdir -p ${ref_dir}/repeat_result/10_rep
mkdir -p ${ref_dir}/repeat_result/15_rep
mkdir -p ${ref_dir}/repeat_result/20_rep
mkdir -p ${ref_dir}/repeat_result/25_rep
mkdir -p ${ref_dir}/repeat_result/50_rep

#genomelist=$(find ${ref_dir}/*repeat.txt | sed 's#.*/##' | awk 'BEGIN{FS="_";OFS="_"}{print $1}')
declare -i x=2
#while [ $x -lt 256 ];
while ((x<=525));
do
    cd ${ref_dir}
    python produce_repeatseq_by_redu_and_txt.py -i ${x}_redundant.fa -l ${x}_repeat.txt -o ${ref_dir}/repeat_result
    cd ${ref_dir}/repeat_result
    mv rep_5.fa ${ref_dir}/repeat_result/5_rep/${x}_rep_5.fa
    mv rep_10.fa ${ref_dir}/repeat_result/10_rep/${x}_rep_10.fa
    mv rep_15.fa ${ref_dir}/repeat_result/15_rep/${x}_rep_15.fa
    mv rep_20.fa ${ref_dir}/repeat_result/20_rep/${x}_rep_20.fa
    mv rep_25.fa ${ref_dir}/repeat_result/25_rep/${x}_rep_25.fa
    mv rep_50.fa ${ref_dir}/repeat_result/50_rep/${x}_rep_50.fa
    x=x+1
done

######
cd ${ref_dir}/repeat_result/5_rep
ls -t *rep_5.fa > 5_rep_list.txt
seqkit stat -j 50 -T --infile-list 5_rep_list.txt > finall_rep5_result.txt
mv finall_rep5_result.txt ../
#######
cd ${ref_dir}/repeat_result/10_rep
ls -t *rep_10.fa > 10_rep_list.txt
seqkit stat -j 50 -T --infile-list 10_rep_list.txt > finall_rep10_result.txt
mv finall_rep10_result.txt ../
######
cd ${ref_dir}/repeat_result/15_rep
ls -t *rep_15.fa > 15_rep_list.txt
seqkit stat -j 50 -T --infile-list 15_rep_list.txt > finall_rep15_result.txt
mv finall_rep15_result.txt ../
######
cd ${ref_dir}/repeat_result/25_rep
ls -t *rep_25.fa > 25_rep_list.txt
seqkit stat -j 50 -T --infile-list 25_rep_list.txt > finall_rep25_result.txt
mv finall_rep25_result.txt ../
#######
cd ${ref_dir}/repeat_result/20_rep
ls -t *rep_20.fa > 20_rep_list.txt
seqkit stat -j 50 -T --infile-list 20_rep_list.txt > finall_rep20_result.txt
mv finall_rep20_result.txt ../
######
cd ${ref_dir}/repeat_result/50_rep
ls -t *rep_50.fa > 50_rep_list.txt
seqkit stat -j 50 -T --infile-list 50_rep_list.txt > finall_rep50_result.txt
mv finall_rep50_result.txt ../

