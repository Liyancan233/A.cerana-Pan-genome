#!/usr/local/bin/bash

#BWA -- burrows-Wheeler Aligner
#GNU core utils
#Java >1.7
#python
#bioawk




mkdir fastq

#ref(assembly) requires fasta format, HiC data requires fastq.gz format

ref=${1}  #ref cat compress
hic_R1=${2} #hic_R1 cat compress
hic_R2=${3}  #hic_R2 cat compress
restriction_enzyme=${4} 



ln -s ${ref} ref_assembly_contigs.fasta
bwa index  ref_assembly_contigs.fasta && echo "*** index is done ***"

cd fastq
ln -s ${hic_R1} hic_R1.fastq.gz 
ln -s ${hic_R2} hic_R2.fastq.gz 

cd ..

python generate_site_positions.py ${restriction_enzyme} ref ref_assembly_contigs.fasta  

awk 'BEGIN{OFS="\t"}{print $1, $NF}' ref_${restriction_enzyme}.txt > genome.chrom.sizes

ln -s {your juicer path}/juicer/CPU/ ./
#ln -s /home/bumblebee/software/juicer/CPU/ ./


bash CPU/juicer.sh -t 50 -z ref_assembly_contigs.fasta -p genome.chrom.sizes -s ${restriction_enzyme} -y ref_${restriction_enzyme}.txt -D CPU && echo "*** juicer is done ***"


#move out the "merged_nodups.txt" file from the obtained "aligned" folder. Then, you could delete all the remaining "aligned" folder.

#change tmp path
export TMPDIR=./


bash ~/software/3d-dna/run-asm-pipeline.sh ref_assembly_contigs.fasta ./aligned/merged_nodups.txt > 3ddna.log 

echo "*** 3d-dna is done ***"

