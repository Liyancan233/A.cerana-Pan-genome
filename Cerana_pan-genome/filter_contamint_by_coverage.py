# -*- coding: utf-8 -*-
import sys
from Bio.Seq import Seq
from Bio import SeqIO
import argparse


parser = argparse.ArgumentParser(description = '\nThis script can remove contaminated sequences from the assembly based on the alignment results，and is generally used with the sequence alignment results', add_help = False, usage = '\npython3 seq_select.py -i [input.fasta] -o [output.fasta] -l [long of contigs ID] -e [contaminated long] -r [coverage rate]')
required = parser.add_argument_group('require parameter')
optional = parser.add_argument_group('optional parameter')
required.add_argument('-i', '--input', metavar = '[input.fasta]', help = 'input，fasta format', required = True)
required.add_argument('-o', '--output', metavar = '[output.fasta]', help = 'output，fasta format', required = True)
required.add_argument('-l', '--long', metavar = '[ID_long.txt]', help = 'he long of contigs ID', required = True)
required.add_argument('-c', '--contaminated', metavar = '[contaminated_list.txt]', help = 'contaminated long', required = True)
required.add_argument('-r', '--rate', metavar = '[number]', help = 'coverage rate', required = True)
optional.add_argument('-h', '--help', action = 'help', help = 'help information')
args = parser.parse_args()
#x=args.rate
'''
list_1=sys.argv[1]#the long of contigs ID
list_2=sys.argv[2]#contaminated long
infile=sys.argv[3]#should filteried fasta
outfile=sys.argv[4]#output file name
rate=sys.argv[5]#rate int
'''
Long_ID=[]
with open(args.long) as f:
	lines1=f.readlines()
	for line1 in lines1:
		listline1=line1.split()
		with open(args.contaminated) as x:
			lines2=x.readlines()
			for line2 in lines2:
				listline2=line2.split()
				if listline1[0] == listline2[0]:
					if int(listline2[1])/int(listline1[1]) > float(args.rate):#change the rate of filteried
						Long_ID.append(listline1[0])
print(Long_ID)
ids = []
for rec in SeqIO.parse(args.input,"fasta"): 
	if rec.id not in Long_ID: 
		ids.append(rec.id)
newfa2=[]
for rec in SeqIO.parse(args.input,"fasta"): 
	if rec.id in ids:
		newfa2.append(rec)
SeqIO.write(newfa2,args.output,"fasta")











