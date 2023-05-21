# -*- coding: utf-8 -*-

import argparse
import sys
from Bio.Seq import Seq
from Bio import SeqIO


parser = argparse.ArgumentParser(description = '\n This script were used to obtain the target gene name', add_help = False, usage = '\npython3 seq_select.py -i [input.fasta] -o [output"name".txt] \n python3 seq_select.py --input [input.fasta]')
required = parser.add_argument_group('required')
optional = parser.add_argument_group('optionial')
required.add_argument('-i', '--input', metavar = '[input.fasta]', help = 'input files，fasta 格式', required = True)
required.add_argument('-o', '--output', metavar = '[output.txt]', help = 'output files，fasta 格式', required = True)

optional.add_argument('-h', '--help', action = 'help', help = 'Information about Help')
args = parser.parse_args()

rep = []
for record in SeqIO.parse(args.input, "fasta"):
    if record.id not in rep:
        rep.append(record.id)
name = args.output
f_write = open(name, 'w')

for rec in rep:
    f_write.write(rec + '\n')
f_write.close()

