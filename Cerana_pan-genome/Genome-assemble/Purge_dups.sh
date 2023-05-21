#! user/loci/bin/bash

asm=${1}
reads=${2}

minimap2 -t 120 -x map-pb $asm $reads | gzip -c - > pb_aln.paf.gz

pbcstat pb_aln.paf.gz
calcuts PB.stat > cutoffs 2> calcults.log
# Split an assembly
split_fa $asm > asm.split
# do a self-self alignment
minimap2 -t 120 -xasm5 -DP asm.split asm.split | gzip -c > asm.split.self.paf.gz
# purge haplotigs and overlap
purge_dups -2 -T cutoffs -c PB.base.cov asm.split.self.paf.gz > dups.bed 2> purge_dups.log
# get the purged primary and haplotigs sequences from draft assembly
get_seqs dups.bed $asm

