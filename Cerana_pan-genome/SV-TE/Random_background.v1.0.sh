#conda base

TE_data=${1}
Observe_data=${2}
Chrom_leng=${3}
bootstrap_NUM=${4}

#unique TE name #HuBei_LG_01     9404    11757   repeat_region   Classification=Unknown

for i in $( seq 1 $bootstrap_NUM ); do echo $i >> tmp.bootstrap_NUM.log;done

cat ${TE_data} | sed  's/Classification=//g' | awk 'BEGIN{FS="\t";OFS="\t"}{print $5}' | sort -n | uniq | sed 1d > TE_name_TMP.txt

cp TE_name_TMP.txt Final_result.txt

#Practical TE count

bedtools intersect -wa -wb -a ${Observe_data} -b ${TE_data} | sed  's/Classification=//g' | awk 'BEGIN{FS="\t";OFS="\t"}{print $11}' > Observe_result.TMP.txt

python Random_background.py -l TE_name_TMP.txt -i Observe_result.TMP.txt > Observe_resulted.TMP.txt

paste Final_result.txt Observe_resulted.TMP.txt > Final_Observe_result.txt

# Random TE count

Random() {
    local i=$1
	bedtools shuffle -i ${Observe_data} -g ${Chrom_leng} -chrom -noOverlapping  > Genomic_background_model-No${i}.TMP.bed 

    bedtools intersect -wa -wb -a Genomic_background_model-No${i}.TMP.bed -b ${TE_data} | sed  's/Classification=//g' | awk 'BEGIN{FS="\t";OFS="\t"}{print $11}' > bootstrap${i}.TMP.txt

    python Random_background.py -l TE_name_TMP.txt -i bootstrap${i}.TMP.txt > Random_result_bootstrap${i}.TMP.txt
}
export -f Random

cat tmp.bootstrap_NUM.log | parallel -j 100 Random

for i in $( seq 1 $bootstrap_NUM ); 
do 
#bedtools shuffle -i ${Observe_data} -g ${Chrom_leng} -chrom -noOverlapping  > Genomic_background_model-No${i}.TMP.bed 

#bedtools intersect -wa -wb -a Genomic_background_model-No${i}.TMP.bed -b ${TE_data} | sed  's/Classification=//g' | awk 'BEGIN{FS="\t";OFS="\t"}{print $11}' > bootstrap${i}.TMP.txt

#python /mnt/cheng/software/lyc_sofeware/SV-TE/Random_background.py -l TE_name_TMP.txt -i bootstrap${i}.TMP.txt > Random_result_bootstrap${i}.TMP.txt

paste Final_result.txt Random_result_bootstrap${i}.TMP.txt > Final_result.TMP.txt

mv Final_result.TMP.txt Final_result.txt

done

rm *TMP*
rm tmp.bootstrap_NUM.log





