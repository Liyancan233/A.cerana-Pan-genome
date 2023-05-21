#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# this script input gff file that maker output ,can output four form file [ exon_stat.txt/exon_size.txt/intron_size.txt/exon_number.txt ] 
#usr gaid [ python this-script input.gff ]
import os
import sys
from itertools import islice
infile_1 = sys.argv[1]
f_write = open("./exon_stat.txt", 'w')
with open (infile_1) as f:
    for line in islice(f, 1, None):
        lin = line.strip().split('\t')
        if lin[0] != '###':
            name = lin[8].split(';')[1].split('=')[-1]
            if lin[2] == 'gene':
                #print '\n',
                #print name, lin[0], lin[3], lin[4],lin[6],
                f_write.write('\n' +  name + '\t' + lin[0] + '\t' + lin[3] + '\t' + lin[4] + '\t' + lin[6]) 
            if lin[2] == 'exon':
                #print lin[3], lin[4],
                f_write.write('\t' + lin[3] + '\t' + lin[4])
f.close()
f_write.close()
##The above program counts the position of exon start and end for each gene
a_write = open("./exon_size.txt", 'w')
with open('exon_stat.txt', 'r') as f:
    for line in f:
        lin = line.strip().split('\t')
        a = len(lin)
        for i in range(6, a, 2):
            exon = int(lin[i]) - int(lin[i-1]) + 1
            #print lin[0], exon
            a_write.write(lin[0] + '\t' + str(exon) + '\n')
f.close()
a_write.close()
#The above program calculates the size of each exon
intron_write = open("./intron_size.txt", 'w')
with open('exon_stat.txt', 'r') as f:
    for line in f:
        lin = line.strip().split('\t')
        a = len(lin)
        if a == 7:
            #print lin[0], '0'
            intron_write.write(lin[0] + '\t' + '0' + '\n')
        if a > 7:
            if lin[4] == '+':
                for i in range(7, a, 2):
                    intron = abs(int(lin[i]) - int(lin[i-1]) - 1)
                    #print lin[0], intron
                    intron_write.write(lin[0] + '\t' + str(intron) + '\n')
            if lin[4] == '-':
                for i in range(8, a, 2):
                    intron = abs(int(lin[i]) + 1 - int(lin[i - 3]))
                    #print lin[0], intron
                    intron_write.write(lin[0] + '\t' + str(intron) + '\n')
f.close()
intron_write.close()
#Calculate the size of each intron
exon_count_write = open("./exon_number.txt", 'w')
with open('exon_stat.txt', 'r') as f:
    for line in islice(f, 1, None):
        lin = line.strip().split('\t')
        a = len(lin)
        n = (a - 5)/2
        #print lin[0], n
        exon_count_write.write(lin[0] + '\t' + str(n) + '\n')
f.close()
exon_count_write.close()









