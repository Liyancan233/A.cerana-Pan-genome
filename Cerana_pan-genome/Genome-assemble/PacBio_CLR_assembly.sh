#!/usr/local/bin/bash

# If the following software is not in the current environment，then add the temporary environment (optional)
export PATH=/{The path to where your software is located}/minimap2-2.13_x64-linux/minimap2:$PATH
export PATH=/{The path to where your software is located}/racon:$PATH
export PATH=/{The path to where your software is located}/bwa-0.7.12/bwa:$PATH

#PB clr data input
PB_file=${1}
#Resequencing data data input
input_file1=${2}
input_file2=${3}
#input raw reads
mkdir canu_racon
mkdir falcon_run
mkdir falcon_racon


echo "*******canu stage*******"
#Correction and canu assembly of raw data

echo "******* start canu *******"
{The path to where your software is located}/canu -d canu_run/ -p Caun maxThreads=100 genomeSize=0.25g minReadLength=2000 minOverlapLength=500 maxMemory=700 corOutCoverage=128 useGrid=false -pacbio-raw  ${PB_file} && echo "**canu is done**"

#The result will appear in /canu_run folder

echo "******* canu done *******"
echo "*******canu start racon*******"

#2.falcon assemble
gzip -d canu_run/Caun.correctedReads.fasta.gz
#Decompress canu's correct fasta.gz
#Configure falcon's configuration file and copy this file to the falcon_run folder
cp fc_bee3.cfg input.fofn falcon_run
cd falcon_run
echo "******* start falcon *******"
fc_run fc_bee3.cfg && echo "** falcon assembly is done ** "
rm -rf 0-rawreads
rm -rf 1-preads_ovl
#assembly result in 2-asm-falcon


cd ../canu_racon
#for falcon dates

#first polish
minimap2 -d ref1.mmi  ../canu_run/Caun.contigs.fasta 
minimap2 -ax map-pb -t 60 ref1.mmi ../canu_run/Caun.correctedReads.fasta > map1.sam  
racon -t 60 ../canu_run/Caun.correctedReads.fasta map1.sam ../canu_run/Caun.contigs.fasta > prefix1.fa && echo "**racon first is done**"
rm map1.sam
# second polish
minimap2 -d ref2.mmi prefix1.fa 
minimap2 -ax map-pb -t 60 ref2.mmi ../canu_run/Caun.correctedReads.fasta > map2.sam 
racon -t 60 ../canu_run/Caun.correctedReads.fasta map2.sam  prefix1.fa  > prefix2.fa && echo "**racon second is done**"
rm map2.sam
# third polish
minimap2 -d ref3.mmi prefix2.fa 
minimap2 -ax map-pb -t 60 ref3.mmi ../canu_run/Caun.correctedReads.fasta > map3.sam 
racon -t 60 ../canu_run/Caun.correctedReads.fasta map3.sam  prefix2.fa  > canu_racon.fa && echo "**racon third is done**"

echo "*******canu racon done *******"
#redundans
echo "******* canu start redundans *******"


redundans.py -v -f canu_racon.fa -o redundans_canu_assembly -t 60 --identity 0.9 --overlap 0.8 --log redundans_canu_assembly.log --noscaffolding --nogapclosing && echo "**redundance is done**"


#canu_resequencing polish
echo "******* canu start pilon *******"
cd redundans_canu_assembly


bwa index scaffolds.reduced.fa


bwa mem -t 60 scaffolds.reduced.fa ${input_file1} ${input_file2} | samtools view -bS - > pctg.bam && echo "** bwa is done**"


samtools sort -f pctg.bam pctg.sorted.bam

#Create the pilon output folder and put the bam from the previous step and the assembled genome into it

mkdir pilonout
cp pctg.sorted.bam scaffolds.reduced.fa pilonout
cd pilonout
samtools index pctg.sorted.bam


java -Xmx600G -jar {The path to where your software is located}/pilon-1.23.jar --genome scaffolds.reduced.fa --bam pctg.sorted.bam --output XX_canu_pilon.fasta --outdir ./ --changes  --fix bases --mindepth 6 --threads 60 echo "**pilon is done**"


echo "******* canu pilon done *******"



echo "******* falcon start racon*******"

cd ../../falcon_racon

#for falcon dates

#first polish
minimap2 -d ref1.mmi ../falcon_run/2-asm-falcon/p_ctg.fa 
minimap2 -ax map-pb -t 60 ref1.mmi ../canu_run/Caun.correctedReads.fasta > map1.sam  
racon -t 60 ../canu_run/Caun.correctedReads.fasta map1.sam ../falcon_run/2-asm-falcon/p_ctg.fa > prefix1.fa && echo "**racon first is done**"
rm map1.sam
# second polish
minimap2 -d ref2.mmi prefix1.fa 
minimap2 -ax map-pb -t 60 ref2.mmi ../canu_run/Caun.correctedReads.fasta > map2.sam 
racon -t 60 ../canu_run/Caun.correctedReads.fasta map2.sam  prefix1.fa  > prefix2.fa && echo "**racon second is done**"
rm map2.sam
# third polish
minimap2 -d ref3.mmi prefix2.fa 
minimap2 -ax map-pb -t 60 ref3.mmi ../canu_run/Caun.correctedReads.fasta > map3.sam 
racon -t 60 ../canu_run/Caun.correctedReads.fasta map3.sam  prefix2.fa  > falcon_racon.fa && echo "**racon third is done**"

echo "******* falcon racon done *******"
#redundans
echo "******* falcon start redundans *******"


redundans.py -v -f falcon_racon.fa -o redundans_falcon_assembly -t 60 --identity 0.9 --overlap 0.8 --log redundans_falcon_assembly.log --noscaffolding --nogapclosing && echo "**redundance is done**"



# resequencing data polish
echo "******* falcon start polish *******"
cd redundans_falcon_assembly


bwa index scaffolds.reduced.fa


bwa mem -t 60 scaffolds.reduced.fa ${input_file1} ${input_file2} | samtools view -bS - > pctg.bam && echo "**bwa is done**"


#samtools view -b -o p_ctg.bam p_ctg.sam
#rm curated.fasta.sam


samtools sort -f pctg.bam pctg.sorted.bam 


mkdir pilonout
cp pctg.sorted.bam scaffolds.reduced.fa pilonout
cd pilonout
samtools index pctg.sorted.bam


java -Xmx600G -jar pilon-1.23.jar --genome scaffolds.reduced.fa --bam pctg.sorted.bam --output XX_falcon_pilon.fasta --outdir ./ --changes  --fix bases --mindepth 6 --threads 60 && echo "**pilon is done**"



echo "******* falcon polish done *******"






















