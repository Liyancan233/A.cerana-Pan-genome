#!/usr/local/bin/bash
#You can run the script by chromosome
#nohup sh combine_filter_indel_snp_vcf.sh /mnt/cheng/lyc/sequence_Pan-genome/ref/lyc_gao_hube_clean_assembly.FINAL.fasta HuBei_LG_02.vcf.gz > log2
#Combined result vcf files
#bcftools concat HuBei_Scf_19_all.raw.gatk.vcf.gz HuBei_Scf_22_all.raw.gatk.vcf.gz -O z -o all.fq.gz


####################################################
# Create a directory, and prepare some path variables for subsequent calls
####################################################
gvcfdir=${1} # folder for g.v.gz
REF=${2} #reference

cd ${gvcfdir}
find ${gvcfdir}/ -name "*.g.vcf.gz" > input.list
mkdir TMP
mv input.list TMP/
cd ${gvcfdir}/TMP
ln -s ${REF} ref.fa
samtools faidx ref.fa
java -Xms2g -Xmx11g -XX:ParallelGCThreads=50 -jar picard.jar CreateSequenceDictionary R=ref.fa O=ref.dict TMP_DIR=./tmp



########################合并GVCF###########################

gatk --java-options "-Xmx80g" CombineGVCFs -R ref.fa \
    --variant input.list \
  -O all.g.vcf.gz --tmp-dir ./

#Convert gvcf to VCF
gatk --java-options "-Xmx80g" GenotypeGVCFs  -R ref.fa \
  -V all.g.vcf.gz \
  -O all.raw.vcf.gz --tmp-dir ./


#############################################################
#Variation result quality control filter, remove low quality variation result
#############################################################

#step1：gatk VariantFiltration(It's just labeled, it's not filtered out) 

gatk --java-options "-Xmx60g" VariantFiltration -R ref.fa -V ./all.raw.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 "  \
    --cluster-window-size 5 --cluster-size 2 --filter-name my_snp_filter -O all.raw.gatk.vcf.gz

#Screen out the mutations that pass the filter
zcat all.raw.gatk.vcf.gz |awk '$0~/#/ || ($7 =="PASS"){print $0}' | gzip - > all.raw.gatked.vcf.gz

#step2：vcfutils.pl
# -w INT    SNP within INT bp around a gap to be filtered [3]
# -W INT    window size for filtering adjacent gaps [10]

vcfutils.pl varFilter -w 5 -W 10 "gzip -dc all.raw.gatked.vcf.gz|" | gzip - > all.varFilter.vcf.gz

#step3：vcftools
#--max-missing Exclude sites on the basis of the proportion of missing data 
#(defined to be between 0 and 1, where 0 allows sites that are completely missing 
#and 1 indicates no missing data allowed).

vcftools --gzvcf all.varFilter.vcf.gz --recode --recode-INFO-all --stdout \
    --maf 0.05  --max-missing 1  --minDP 4  --maxDP 1000  \
    --minQ 30 --minGQ 80 --min-alleles 2  --max-alleles 2 | gzip - > all.clean.vcf.gz

#分离indel与SNP 到不同的文件
vcftools --remove-indels --recode --recode-INFO-all --gzvcf all.clean.vcf.gz --stdout |gzip - > all.clean.snp.vcf.gz
vcftools --keep-only-indels  --recode --recode-INFO-all --gzvcf all.clean.vcf.gz --stdout |gzip - > all.clean.indel.vcf.gz
