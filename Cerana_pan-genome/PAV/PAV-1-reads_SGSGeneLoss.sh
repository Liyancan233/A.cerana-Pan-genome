#!/usr/local/bin/bash

ref_dir=${1}  #For example/mnt/cheng/lyc/genome/ref.fa
reads_dir=${2} #For example /mnt/cheng/lyc/genome
out_put_dir=${3} #For example /mnt/cheng/lyc/genome
#requrie_dir=${4} 
gff_dir=${4}

genomelist=$(find ${reads_dir}/*_1.fastq.gz | sed 's#.*/##' | awk 'BEGIN{FS="_";OFS="_"}{print $1,$2,$3}')
#fastq name is XX_YY_ZZ_1(2).fastq.gz (for example FuJian_LongYan_357357_1.fastq.gz)

#mkdir -p 6_place_contigs
cd ${out_put_dir}
mkdir result
cp $ref_dir ./ref.fa
bwa index -p ref.fa ${ref_dir}
samtools faidx ref.fa


for x in $genomelist;
do
        loci=$(echo $x | awk -v FS="_" '{print $1}') #fastq name is XX_YY_ZZ_1(2).fastq.gz (for example FuJian_LongYan_357357_1.fastq.gz)
        sample=$(echo $x | awk 'BEGIN{FS="_";OFS="_"}{print $2,$3}')
        Id=$(echo $x | awk -v FS="_" '{print $3}')
        cd ${out_put_dir}
        mkdir -p ${out_put_dir}/TMP/$x
        cd ${out_put_dir}/TMP/$x

        ln -s ${ref_dir} ref.fa

        /mnt/cheng/software/bwa-0.7.12/bwa mem ref.fa -t 40 -R "@RG\tID:${Id}\tLB:${loci}\tSM:${sample}\tPL:ILLUMINA" ${reads_dir}/${x}_1.fastq.gz ${reads_dir}/${x}_2.fastq.gz | /mnt/cheng/software/samtools-1.7/samtools view -@ 40 -bS -f 2 - | /mnt/cheng/software/samtools-1.7/samtools sort -@ 40 > test_read.sorted.bam
        #/mnt/cheng/software/samtools-1.7/samtools sort -@ 40 -o test_read.sorted.bam test_read.bam
        /mnt/cheng/software/samtools-1.7/samtools index -@ 40 test_read.sorted.bam
        mkdir -p ${out_put_dir}/result/$x
        java -Xmx80g -jar /mnt/cheng/software/new_lyc_software/SGSGeneLossv0.1/SGSGeneLoss.v0.1.jar bamPath=${out_put_dir}/TMP/${x}/ bamFileList=test_read.sorted.bam gffFile=${gff_dir} outDirPath=${out_put_dir}/result/${x}/ chromosomeList=all minCov=2 lostCutoff=0.2
done        
