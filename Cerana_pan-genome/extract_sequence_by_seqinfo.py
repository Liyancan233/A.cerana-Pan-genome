# -*- coding: utf-8 -*-

import argparse
import sys
from Bio.Seq import Seq
from Bio import SeqIO


parser = argparse.ArgumentParser(description = '\nThis script is used to intercept sequences at specific locations in the genome, requiring additional input of a list file recording the intercepted sequence information', add_help = False, usage = '\npython3 seq_select.py -i [input.fasta] -o [output.fasta] -l [list]\npython3 seq_select.py --input [input.fasta] --output [output.fasta] --list [list]')
required = parser.add_argument_group('Required parameters')
optional = parser.add_argument_group('Optional parameters')
required.add_argument('-i', '--input', metavar = '[input.fasta]', help = 'input files，fasta format', required = True)
required.add_argument('-o', '--output', metavar = '[output.fasta]', help = 'output files，fasta format', required = True)
required.add_argument('-l', '--list', metavar = '[list]', help = 'Record "sequence ID/ sequence start position/sequence end position/separated by tab', required = True)
optional.add_argument('-h', '--help', action = 'help', help = 'information')
args = parser.parse_args()

newfa2=[]       
dic={}
O=''

with open(args.list) as f:
		lines=f.readlines() #The output will have a \n after each line
		for line in lines:
				list1=[]
				listline=line.strip().split('\t')  #Splitting the line is equivalent to removing the final \n
				if int(listline[1]) > int(listline[2]):
						O = listline[2]
						listline[2] = listline[1]
						listline[1] = O
				if listline[0] in dic.keys():
						dic[listline[0]].append(tuple(listline[1:3]))
				else:
						list1.append(tuple(listline[1:3]))
						dic[listline[0]]=list1


print(dic) 

for record in SeqIO.parse(args.input, "fasta"):
	
	if record.id in dic.keys():
		for i in dic[record.id]:     
			a=record.seq[0:int(int(i[0])-1)] + "K"*int(int(i[1])-int(i[0])+1) + record.seq[int(i[1]):]
			#a=record.seq(int(num1)-1:int(num2):)
		#b=str(a.seq).replace("K",'')
			record.seq=a
		record.seq=Seq(str(record.seq).replace("K",'N'))
		newfa2.append(record)
	else:
		newfa2.append(record)

SeqIO.write(newfa2, args.output, "fasta")

