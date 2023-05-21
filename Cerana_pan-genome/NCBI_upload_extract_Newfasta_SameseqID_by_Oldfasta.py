# -*- coding: utf-8 -*-


import argparse
import sys
from Bio.Seq import Seq
from Bio import SeqIO


parser = argparse.ArgumentParser(description = '\n This script is used for NCBI to upload contaminated fragments returned from genomic data.', add_help = False, usage = '\npython3 seq_select.py -N [New.fasta] -O [Old.fasta] -out [Newfasta_with_Oldseqname.fasta] \n python3 seq_select.py --input [input.fasta]')
required = parser.add_argument_group('Required parameters')
optional = parser.add_argument_group('Optional parameters')
required.add_argument('-N', '--newfatsa', metavar = '[newfasta.fasta]', help = 'input files，fasta format', required = True)
required.add_argument('-O', '--oldfasta', metavar = '[oldfasta.fasta]', help = 'output files，fasta format', required = True)
required.add_argument('-out', '--output', metavar = '[output.fasta]', help = 'output files，fasta 格式', required = True)

optional.add_argument('-h', '--help', action = 'help', help = '帮助信息')
args = parser.parse_args()
#NCBI auto clean fasta result Name forme lcl|HuBei_LG_01 Apis cerana
rep = []
for record in SeqIO.parse(args.oldfasta, "fasta"):
    if record.id not in rep:
        rep.append(record.id)
new_seq = []
for record in SeqIO.parse(args.newfatsa, "fasta"):
    new_seqname = record.id.strip().split('|')[1]
    if new_seqname in rep:
        print (new_seqname)
        new_seq.append(record)
SeqIO.write(new_seq, args.output, "fasta")

