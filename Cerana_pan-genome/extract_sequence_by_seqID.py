# -*- coding: utf-8 -*-
import sys
from Bio.Seq import Seq
from Bio import SeqIO

infile1=sys.argv[1]
infile2=sys.argv[2]
outfile=sys.argv[3]
#use this_py input1(fasta) input2(the ID that we wang to filter) output(remain fasta)

id = []
for line in open(infile2):
    id.append(line.rstrip().strip('"'))

newfa2=[]
for rec in SeqIO.parse(infile1,"fasta"):
    if rec.id in id:
        newfa2.append(rec)
SeqIO.write(newfa2, outfile, "fasta")


