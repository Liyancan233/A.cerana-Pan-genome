#!/usr/local/bin/bash
ref_dir=${1}
rm *clean.fasta
rm redundant.fa ref.fa list*
ls -t lyc_* > nonredu.txt
ls -t *redundant.fa > redu.txt
seqkit stat -j 50 -T --infile-list nonredu.txt > nonredu_result.txt
seqkit stat -j 50 -T --infile-list redu.txt > redu_result.txt
paste nonredu_result.txt redu_result.txt | sed '1d' | awk 'BEGIN{FS="\t";OFS="\t"}{print $9,$4,$5,$7,$4-$12,$5-$13,($5-$13)/($4-$12),$12,$13,$15}' > nonredu_unique_common.txt
awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$3,$6,$9,$2,$5,$8,$4,$7,$10}' nonredu_unique_common.txt > size_number_leng.txt
mkdir -p result
mv nonredu_unique_common.txt result/
mv size_number_leng.txt result/
rm nonredu.txt redu.txt nonredu_result.txt redu_result.txt
#sh produce_repeatseq_by_redu_and_txt.sh ${ref_dir}



