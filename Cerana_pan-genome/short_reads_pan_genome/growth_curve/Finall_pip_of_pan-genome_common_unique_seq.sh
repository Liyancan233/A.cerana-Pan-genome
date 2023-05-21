#!/usr/local/bin/bash
# Run the script in the 5_finally_contigs folder generated in the previous step, using the absolute path for the relevant scripts used in this run
ref_dir=${1}  #for rxample /mnt/cheng/lyc/genome/ref.fa, Take care to create an individual sequence as initial "ref.fa" when using this script

genomelist=$(find ${ref_dir}/*clean.fasta | sed 's#.*/##' | awk 'BEGIN{FS="_";OFS="_"}{print $0}' | sort --random-sort)
y=1
touch redundant.fa
touch repeat.txt
cp ref.fa root_ref.fa
#BLASTN
#/mnt/cheng/software/ncbi-blast-2.11.0+/bin/makeblastdb -in ${ref_dir} -dbtype nucl -out rf_Id#
#/mnt/cheng/software/ncbi-blast-2.11.0+/bin/blastn -task blastn -db ref_Id -query ${ref_dir} -evalue 1e-5 -best_hit_overhang 0.25 -perc_identity 0.8 -max_hsps 1 -num_threads 30 -outfmt '6 qseqid sseqid pident slen qlen length qstart qend sstart send mismatch gapopen gaps evalue bitscore' -out output_blastnt.txt
#nohup python /mnt/cheng/software/lyc_sofeware/blastn.py --alignment_path output_blastnt.txt --input ${ref_dir} --identity 90 --coverage 0.90 --output blastn.fasta &#
for x in $genomelist;
do
    y=$[$y+1];
    nucmer --maxmatch -l 31 -c 100 -t 80 -p ${y}_deltaFile ref.fa ${x}
    show-coords -H -T -l -c -o ${y}_deltaFile.delta > ${y}_coordsFile
    python produce_rep_nonredu_and_commonSeq_3.0.py --alignment_path ${y}_coordsFile --input_ref ref.fa --input_qry ${x} --input_oldredundant redundant.fa --identity 0.9 --coverage 0.80 --output1 ref.fa --output_newredundant redundant.fa
    cp redundant.fa ${y}_redundant.fa
    cp ref.fa lyc_${y}_ref.fasta
    rm ${y}_deltaFile.delta ${y}_coordsFile
    cp new_repeat.txt repeat.txt
    mv new_repeat.txt ${y}_repeat.txt
done

echo " step 1 is done "

#********************************************************************************************

# produce un-unique un-reference sequence with repeats 5,10,15,20,25
sh produce_repeatseq_by_redu_and_txt.sh $ref_dir

echo " step 2 is done "
#**********************************************************************************************

#produce result

sh result_nonredu_common_statelist.sh ${ref_dir}


echo " step 3 is done "

#The repeated results are in the result/ folder. It is recommended to use this script several times to evaluate the pangenomic growth trend (100 repeats are recommended).


