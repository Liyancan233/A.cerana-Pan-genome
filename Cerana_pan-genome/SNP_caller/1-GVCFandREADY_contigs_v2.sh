#!/usr/local/bin/bash

ref_dir=${1}  #e.g./mnt/cheng/lyc/genome/ref.fa
genome_dir=${2} #e.g./mnt/cheng/lyc/genome
#name=${3} #e.g.hube
out_put_dir=${3} #e.g./mnt/cheng/lyc/genome

genomelist=$(find ${genome_dir}/*_1.fastq.gz | sed 's#.*/##' | awk 'BEGIN{FS="_";OFS="_"}{print $1,$2,$3}')


mkdir -p ${out_put_dir}/GVCF
#mkdir -p ${out_put_dir}/Ready_contigs
#mkdir -p ${out_put_dir}/alignedmate
for x in $genomelist;
do
	loci=$(echo $x | awk -v FS="_" '{print $1}')
	sample=$(echo $x | awk 'BEGIN{FS="_";OFS="_"}{print $2,$3}')
        Id=$(echo $x | awk -v FS="_" '{print $3}')
	cd ${genome_dir}
	java -Xms2g -Xmx11g -jar /home/liyancan/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/trimmomatic.jar PE ${x}_1.fastq.gz ${x}_2.fastq.gz -baseout ${x}.fq.gz ILLUMINACLIP:/home/liyancan/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/TruSeq2-PE.fa:2:30:10:2:keepBothReads -threads 50 LEADING:3 MINLEN:70
	mkdir -p ${genome_dir}/TMP/$x
	#mkdir -p ${out_put_dir}/alignedmate/$x
	cd ${genome_dir}/TMP/$x
	ln -s ${ref_dir} ref.fa
	/mnt/cheng/software/bwa-0.7.12/bwa index ref.fa
	#bwa aln ref.fa ${genome_dir}/${name}_${x}_1P.fq.gz -t 80 > ${name}_${x}_1.fq.sai
	#bwa aln ref.fa ${genome_dir}/${name}_${x}_2P.fq.gz -t 80 > ${name}_${x}_2.fq.sai 
	#bwa sampe ref.fa -r "@RG\tID:<ID>\tLB:<${name}>\tSM:<${x}>\tPL:ILLUMINA" ${name}_${x}_1.fq.sai ${name}_${x}_2.fq.sai ${genome_dir}/${name}_${x}_1P.fq.gz ${genome_dir}/${name}_${x}_2P.fq.gz > test_read.sam
	samtools faidx ref.fa
	/mnt/cheng/software/bwa-0.7.12/bwa mem ref.fa -t 40 -R "@RG\tID:${Id}\tLB:${loci}\tSM:${sample}\tPL:ILLUMINA" ${genome_dir}/${x}_1P.fq.gz ${genome_dir}/${x}_2P.fq.gz | samtools view -S -b - > test_read.bam
	/home/liyancan/miniconda3/bin/java -Xms2g -Xmx11g -XX:ParallelGCThreads=50 -jar /home/liyancan/miniconda3/share/picard-2.25.7-0/picard.jar CreateSequenceDictionary R=ref.fa O=ref.dict TMP_DIR=./tmp
	#picard ReorderSam I= test_read.sam O= test_reads.reordered.sam SD=ref.fa TMP_DIR=./tmp
	#samtools view -bS test_reads.reordered.sam -o test.bam
	/home/liyancan/miniconda3/bin/java -Xms2g -Xmx11g -XX:ParallelGCThreads=50 -jar /home/liyancan/miniconda3/share/picard-2.25.7-0/picard.jar SortSam INPUT=test_read.bam OUTPUT=sorted.bam SORT_ORDER=coordinate TMP_DIR=./tmp
	#time samtools sort -@ 40 -m 11G test_read.bam -o sorted.bam
	/home/liyancan/miniconda3/bin/java -Xms2g -Xmx11g -XX:ParallelGCThreads=50 -jar /home/liyancan/miniconda3/share/picard-2.25.7-0/picard.jar MarkDuplicates REMOVE_DUPLICATES= false MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 INPUT= sorted.bam OUTPUT= sorted_dedup.bam METRICS_FILE= sorted_dedup.metrics TMP_DIR=./tmp
	
	samtools index sorted_dedup.bam
	echo"**succeed!!!**"
	nohup gatk --java-options '-Xmx20G' HaplotypeCaller -R ref.fa -O ${x}_output.g.vcf.gz -I sorted_dedup.bam -ERC GVCF & 
	#cp ${name}_${x}_output.g.vcf.gz ${out_put_dir}/GVCF
	rm test_read.bam sorted.bam
	cd ${genome_dir}
	rm ${x}_1P.fq.gz ${x}_1U.fq.gz ${x}_2P.fq.gz ${x}_2U.fq.gz
done


####
#find ./ -name "*.g.vcf.gz" > input.list
#sh /mnt/cheng/software/lyc_sofeware/index_gvcf.sh ./
#/home/liyancan/miniconda3/bin/java -Xms2g -Xmx11g -XX:ParallelGCThreads=50 -jar /home/liyancan/miniconda3/share/picard-2.25.7-0/picard.jar CreateSequenceDictionary R=ref.fa O=ref.dict TMP_DIR=./tmp
##gatk --java-options '-Xmx500G' CombineGVCFs -R /gatk/input/reference/reference.fa --variant input.list -O combined.g.vcf

#gatk --java-options '-Xmx500G' GenotypeGVCFs  -R /gatk/input/reference/reference.fa -V /gatk/output/gvcf/combined.g.vcf -O /gatk/output/vcf/genome.vcf

#/mnt/cheng/software/gatk-4.1.4.0/gatk GenotypeGVCFs -R ${ref} -V ${data}/bam/${n}.g.raw.vcf -O ${data}/bam/${n}.raw.vcf

#/mnt/cheng/software/gatk-4.1.4.0/gatk SelectVariants -V ${data}/bam/${n}.raw.vcf -select-type SNP -O ${data}/bam/${n}.select.snp.vcf

#/mnt/cheng/software/gatk-4.1.4.0/gatk SelectVariants -V ${data}/bam/${n}.raw.vcf -select-type INDEL -O ${data}/bam/${n}.select.indel.vcf

#/mnt/cheng/software/gatk-4.1.4.0/gatk VariantFiltration -V ${data}/bam/${n}.select.snp.vcf --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "PASS" -O ${data}/bam/${n}.filter.select.snp.vcf
###









