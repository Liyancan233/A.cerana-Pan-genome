#!/usr/local/bin/bash
SGSG_result_file=${1}  #For example /mnt/cheng/lyc/genome/result
geneID_list=${2} #combine_seqID.txt,this file is used in up step
output_folder=${3} 
#mkdir -p ${SGSG_result_file}/temp
samplelist=$(ls ${SGSG_result_file} | sed 's#\:##' | awk 'BEGIN{FS="_";OFS="_"}{print $1,$2,$3}')
chromID=$(cat $geneID_list | awk 'BEGIN{FS="_";OFS="_"}{print $0}')


for x in $samplelist;
do
    A=1
    B=1
    #cd ${SGSG_result_file}/$x
    for y in ${chromID};
    do
        cd ${SGSG_result_file}/$x
        if [ $B == $A ];then
            cat ${y}.excov | awk 'BEGIN{FS=",";OFS="\t"}{print $1,$2,$4,$5,$3}'
        else
            cat ${y}.excov | sed '1d' | awk 'BEGIN{FS=",";OFS="\t"}{print $1,$2,$4,$5,$3}'
        fi
        A=$[$A+1];
    done > ${output_folder}/${x}.tsv
    #cd ${output_folder}
    #paste ${x}.tsv
    echo $x
done > ${output_folder}/samplelist.txt

cd ${output_folder}
touch boost.txt
z=1
b=1 
for x in $samplelist;
do    
    if [ $b == $z ];then
        cat ${x}.tsv | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,$4,$5}' > qry.txt
    else
        cat ${x}.tsv | awk 'BEGIN{FS="\t";OFS="\t"}{print $5}' > qry.txt
    fi
    paste boost.txt qry.txt > pro_boost.txt
    mv pro_boost.txt boost.txt
    z=$[$z+1];
done

# sed 's/PRESENT/1/g' 525.txt > step1_525.txt
# sed 's/LOST/0/g' step1_525.txt > step2_525.txt
