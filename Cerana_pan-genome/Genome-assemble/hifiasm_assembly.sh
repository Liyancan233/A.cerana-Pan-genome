#ccs rawreads data (bam) for three generations of HiFi requires pre-conversion using pcbio's ccs software
#conda
ccs m64064_191201_134422.subreads.0--0.bam hubei.reads.fastq.gz
#hifiasm
# Run on test data (use -f0 for small datasets)

hifiasm -o test -t4 -f0 chr11-2M.fa.gz 2 > test.log
#change to fasta
awk '/^S/{print ">"$2;print $3}' test.p_ctg.gfa > test.p_ctg.fa  # get primary contigs in FASTA
