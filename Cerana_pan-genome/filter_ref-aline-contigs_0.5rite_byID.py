# -*- coding: utf-8 -*-
import sys
from Bio.Seq import Seq
from Bio import SeqIO


list_1=sys.argv[1]#the long of contigs ID
list_2=sys.argv[2]#extract the contigs len
infile=sys.argv[3]#should filteried fasta
outfile=sys.argv[4]#output file name
rate=sys.argv[5]#rate int
Long_ID=[]
with open(list_1) as f:
	lines1=f.readlines()
	for line1 in lines1:
		listline1=line1.split()
		with open(list_2) as x:
			lines2=x.readlines()
			for line2 in lines2:
				listline2=line2.split()
				if listline1[0] == listline2[0]:
					if int(listline2[1])/int(listline1[1]) > float(rate):#change the rate of filteried
						Long_ID.append(listline1[0])
print(Long_ID)
ids = []
for rec in SeqIO.parse(infile,"fasta"): 
	if rec.id not in Long_ID: 
		ids.append(rec.id)
newfa2=[]
for rec in SeqIO.parse(infile,"fasta"): 
	if rec.id in ids:
		newfa2.append(rec)
SeqIO.write(newfa2, outfile, "fasta")

