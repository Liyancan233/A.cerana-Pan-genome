#!/usr/local/bin/bash

ref_dir=${1}  #For example, /mnt/cheng/lyc/genome/ref.fa
ref_mitochondrial=${2}
genome_dir=${3} #For example, /mnt/cheng/lyc/genome
#name=${4} #For examplehube
out_put_dir=${4} #For example, /mnt/cheng/lyc/genome
genomelist=$(find ${genome_dir}/*_1.fastq.gz | sed 's#.*/##' | awk 'BEGIN{FS="_";OFS="_"}{print $1,$2,$3}')
export BLASTDB='/mnt/huang/nt'
cd ${out_put_dir}
mkdir -p 5_finally_contigs
mkdir -p 1_frist_assembly
mkdir -p 2_non-reference_contigs
mkdir -p 3_redundant_and_filter-over-500_contigs
mkdir -p 4_remove_contamint
#mkdir -p 6_place_contigs
cd ${genome_dir}

bwa index -p ref ${ref_dir}
makeblastdb -in ${ref_dir} -dbtype nucl -out ref_alt_Id
makeblastdb -in ${ref_mitochondrial} -dbtype nucl -out ref_mito_Id


for x in $genomelist;
do
	cd ${genome_dir}
	mkdir -p ${genome_dir}/TMP/$x
	cd ${genome_dir}/TMP/$x
	#This step obtains unaligned reads, including double-ended, single-ended fastq format data
	bwa mem -t 80 ${genome_dir}/ref ${genome_dir}/${x}_1.fastq.gz ${genome_dir}/${x}_2.fastq.gz > alignment.sam
	samtools fastq -f 12 -@ 80 alignment.sam -1 R1_Unalignedmate.fq -2 R2_Unalignedmate.fq
	samtools fastq -f 68 -F 8 -@ 80 alignment.sam > R1_alignedmate.fq
	samtools fastq -f 132 -F 8 -@ 80 alignment.sam > R2_alignedmate.fq
	samtools view -b -f 8 -F 4 -@ 80 alignment.sam > alignedmate_ref.bam
	megahit -t 80 -1 R1_Unalignedmate.fq -2 R2_Unalignedmate.fq -r R1_alignedmate.fq, R2_alignedmate.fq -o sample_1
	#Reference genes were aligned using mucmer, similarity over 0.9, and align proportion over 0.8 were excluded
	cd megahit_out
	cp final.contigs.fa ${out_put_dir}/1_frist_assembly/${x}_contigs.fa
	blastn -task blastn -db ${genome_dir}/ref_alt_Id -outfmt 6 -out output_blastnt.txt -max_hsps 1 -num_threads 96 -max_target_seqs 1 -query final.contigs.fa
	perl fasta_size_report.pl -f final.contigs.fa | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2 }' > long.txt #Get the contigs long
	awk 'BEGIN{FS="\t";OFS="\t"}{if($3>90) print $1,$4}' output_blastnt.txt > aling_leng.txt #identity>90. Get the contigs mapped identity
	python lyc_filter_contamint_by_coverage.py -l long.txt -c aling_leng.txt -i final.contigs.fa -o filter_indivi.fasta -r 0.8 #≥50% coverage
	cp filter_indivi.fasta ${out_put_dir}/2_non-reference_contigs/${x}_nonref_indivi.fasta
	#contigs with lengs less than 500bp were filtered and redundances were removed
	perl filter_fasta_by_size.pl -f filter_indivi.fasta -s 500 -u 0.5_long_500.fasta 
	python prepend_to_fasta_header.py -i 0.5_long_500.fasta -o head_out.fasta -p ${x} #add tag at individual
	cd-hit-est -c 0.9 -G 0 -aL 0.90 -AL 500 -aS 0.9 -T 30 -M 0 -i head_out.fasta -o redundant_head_out.fasta #Redundancy removal
	cp redundant_head_out.fasta ${out_put_dir}/3_redundant_and_filter-over-500_contigs/${x}_redundant_head_out.fasta
	#Remove contaminants, archaea, bacteria, human, UniVec_Core, viral
	mkdir contaminant
	cd contaminant
	kraken2 --db super_db/ --report name.output.txt --use-names --threads 60 --unclassified-out kraken2_redundant_head_out.fasta ../redundant_head_out.fasta > classified.txt
	#Use blastn's NT library for deep filtering to remove non-insect genes
	#blastn -task blastn -db nt -evalue 1e-5 -best_hit_overhang 0.25 -perc_identity 0.8 -outfmt '6 qaccver saccver pident length mismatch ssciname qstart qend sstart send evalue bitscore' -out output_blastnt.txt -max_hsps 1 -num_threads 96 -query kraken2_redundant_head_out.fasta
	#perl fasta_size_report.pl -f kraken2_redundant_head_out.fasta | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2 }' > long.txt 
	#awk 'BEGIN{FS="\t";OFS="\t"}{if($3>90) print}' output_blastnt.txt > output_blastnt.txt.goodmatch.txt 
	#perl /mnt/huang/contaminant_filtering/script/get_speciesname_from_blastn.pl output_blastnt.txt.goodmatch.txt #Extract the specified species
	#python /mnt/huang/contaminant_filtering/script/get_species_placement_in_NCBI.py output_blastnt.txt.goodmatch.txt.txt title.txt #Retrieves all species classification paths
	#paste -d /t output_blastnt.txt.goodmatch.txt title.txt > final_blastnt.txt
	#sed -i -e '/Insecta/d' final_blastnt.txt #Leave all non-insect contigs
	#awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$4}' final_blastnt.txt > aling_leng.txt 
	#python lyc_filter_contamint_by_coverage.py -l long.txt -c aling_leng.txt -i kraken2_redundant_head_out.fasta -o ${x}_remove_contaminant.fasta -r 0.5 #≥50% coverage
	cp kraken2_redundant_head_out.fasta ${x}_remove_contaminant.fasta
	cp ${x}_remove_contaminant.fasta ${out_put_dir}/4_remove_contamint
	#The next step is to remove mitochondrial contamination
	mkdir remove_mitochondrial
	cd remove_mitochondrial
	blastn -task blastn -db ${genome_dir}/ref_mito_Id -outfmt 6 -out output_blastnt.txt -max_hsps 1 -num_threads 96 -max_target_seqs 1 -query ../${x}_remove_contaminant.fasta
	perl fasta_size_report.pl -f ../${x}_remove_contaminant.fasta | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2 }' > long.txt 
	awk 'BEGIN{FS="\t";OFS="\t"}{if($3>90) print $1,$4}' output_blastnt.txt > aling_leng.txt
	python lyc_filter_contamint_by_coverage.py -l long.txt -c aling_leng.txt -i ../${x}_remove_contaminant.fasta -o ${x}_clean_indivi.fasta -r 0.8 #≥80% coverage
	cp ${x}_clean_indivi.fasta ${genome_dir}/TMP/$x
	cp ${x}_clean_indivi.fasta ${out_put_dir}/5_finally_contigs
	cd ${genome_dir}/TMP/$x
	rm alignment.sam
	#/mnt/cheng/software/bowtie2-2.2.6/bowtie2-build ${x}_clean_indivi_fasta contig_Id
	#/mnt/cheng/software/bowtie2-2.2.6/bowtie2 -x contig_Id -U R1_alignedmate.fq,R2_alignedmate.fq -S readtocontig.sam
	#samtools view -h -F 2304 readtocontig.sam | samtools sort -n -O bam | /mnt/cheng/software/bedtools2/bin/bedtools bamtobed -i stdin | awk '{OFS="\t"} {print $4,$1,$6,$2,$3}' | sed -e 's/\/[1-2]//g' |sort > readtocontig.txt
	#samtools view -H alignedmate_ref.bam | cat - <"(awk 'NR==FNR{ a[$1]; next }$1 in a{ print $0 ; delete a[$1]; next }' readtocontig.txt <(samtools view alignedmate_ref.bam ))" | samtools sort -n -O bam | /mnt/cheng/software/bedtools2/bin/bedtools bamtobed -i stdin | awk '{OFS="\t"}{print $4,$1,$6,$2,$3}' | sed -e 's/\/[1-2]//g' | sort > pass_mates.txt
	#join -j 1 readtocontig.txt pass_mates.txt > ${x}_mates_region.txt
	#cp ${x}_mates_region.txt ${out_put_dir}/6_place_contigs
done
