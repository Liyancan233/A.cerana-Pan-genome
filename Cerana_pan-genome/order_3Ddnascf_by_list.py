# -*- coding: utf-8 -*-
import sys
from Bio.Seq import Seq
from Bio import SeqIO
import argparse


parser = argparse.ArgumentParser(description = '\nThis script sorts the chromosomes of the specified order according to the list information, which is used for chromosome order recalipation after 3ddna assembly', add_help = False, usage = '\npython order_3Ddnascf_by_list.py -i [input.fasta] -o [output.fasta] -l [order list]')
required = parser.add_argument_group('required parameter')
optional = parser.add_argument_group('optional parameter')
required.add_argument('-i', '--input', metavar = '[input.fasta]', help = 'input，fasta format', required = True)
required.add_argument('-o', '--output', metavar = '[output.fasta]', help = 'output，fasta format', required = True)
required.add_argument('-l', '--list', metavar = '[correct_IDlist.txt]', help = 'correct order of contigs ID', required = True)
# example list : oldname -> newname
# HiC_scaffold_1	HiC_scaffold_6
# HiC_scaffold_2	HiC_scaffold_12
# HiC_scaffold_3	HiC_scaffold_7
optional.add_argument('-h', '--help', action = 'help', help = 'gelp inforamation')
args = parser.parse_args()
ID = []
dic={}

with open(args.list) as f:
    for line in f:
        lin = line.strip().split('\t')
        if lin[0] not in ID:
            ID.append(lin[0])
            dic[lin[0]] = lin[1]
            
newfa2=[]

for rec in SeqIO.parse(args.input,"fasta"):

    if rec.id in ID:
        a = dic[rec.id]
        #print(dic[rec.id])
        rec.id = a
        rec.description = a
        print(rec.id)
        print(rec.description)
        newfa2.append(rec)
    else:
        newfa2.append(rec)
SeqIO.write(newfa2,args.output,"fasta")
