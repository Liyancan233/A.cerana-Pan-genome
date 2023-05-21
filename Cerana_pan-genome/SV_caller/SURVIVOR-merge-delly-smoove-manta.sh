#!/usr/local/bin/bash
#requrie Python2
smoove_dir=${1} # output of Individual smoove call sv [3-genotype_each_sample_at_sites (or) 1-step-file ]
delly_dir=${2} # delly result  [Individual_sample_genotypy]
manta_dir=${3} # manta result  [/mnt/huang/lyc/SV_caller/SV_manta/1-all-sample-result/RESULT]
output_dir=${4} # 
#Merge_sv_sites=${5}

# we select a floder to run for example delly ID name HaiNan_HaiKou_335335_geno.vcf.gz
# step 1 merge diff software output with a same sample

#*********************************************  step one *******************************************************************
mkdir -p ${output_dir}/Step1-merge-diff-software
mkdir -p ${output_dir}/Step2-merge-diff-sample
echo "**********************************************step 1 start***********************************************"
#after filter quality of sites we shoud another remove pop-genotyped sites likes './.','0/0', we only used individual information
genomelist=$(find ${delly_dir}/*.vcf.gz | sed 's#.*/##' | awk 'BEGIN{FS="_";OFS="_"}{print $1,$2,$3}')
y=1
for x in $genomelist;
do
    cd ${output_dir}/Step1-merge-diff-software
    #smoove samplle SR<3 are remove
    bcftools view -e "INFO/SR<3" -Ov ${smoove_dir}/${x}-smoove.genotyped.vcf.gz | grep -v '0/0:' | grep -v '\.\/\.\:' > ${x}-3-smoove-cleand.vcf
    #delly sample LowQual are remove
    bcftools view -f "PASS" -Ov ${delly_dir}/${x}_geno.vcf.gz | grep -v '0/0:' | grep -v '\.\/\.\:' > ${x}-2-delly-cleand.vcf
    #manta sample uncopress
    bcftools view -OV ${manta_dir}/${x}_diploidSV.vcf.gz | grep -v '0/0:' | grep -v '\.\/\.\:' > ${x}-1-manta-cleand.vcf
    #start to merge
    ls ${x}* > sample_files
    #Two adjacent SVs were combined as a single SV if the distance between start coordinate of one SV and end coordinate of the other SV was less than 500 bp.
    /mnt/cheng/software/SURVIVOR-master/Debug/SURVIVOR merge sample_files 500 2 1 1 0 50 ${x}_merged.vcf && echo " $x is done,total number $y "
    #500 means maximum allowed distance of 1kb ,as measured pairwise between breakpoints (begin1 vs begin2, end1 vs end2);
    #2 means ask SURVIVOR only to report calls supported by 2 callers;
    #1 agree on the type;
    #1 agree on the strand;
    #0 this ie invalid parameter
    #50 min SV len
    y=$[$y+1]
    mv *merged.vcf ${output_dir}/Step2-merge-diff-sample
done
echo "*********************************************step 1 done*********************************************"

echo "*********************************************step 2 start***********************************************"

cd ${output_dir}/Step2-merge-diff-sample
ls *merged.vcf >sample_files
/mnt/cheng/software/SURVIVOR-master/Debug/SURVIVOR merge sample_files 500 1 0 0 0 50 SURVIVOR-Final_merged.vcf

echo "***************************************SURVIVOR is done !!!***********************************"


#I used SURVIVOR-Final_merged.vcf site  base on delly data EXCEL manual script!!


#################################################################################################################
