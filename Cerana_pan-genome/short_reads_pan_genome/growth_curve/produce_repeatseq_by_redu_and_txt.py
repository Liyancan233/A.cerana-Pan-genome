# -*- coding: utf-8 -*-

#导入模块，初始传递命令、变量等
import argparse
import sys
import os
from Bio.Seq import Seq
from Bio import SeqIO

parser = argparse.ArgumentParser(description = '\n该脚本用于将输入的序列按序列名称的第一个名字（HuBei）为生成的文件名称（>HuBei_ShenNongJ_SRR10549768k119_2959，且分隔符为"_"）将序列分割,并输出到result文件夹内', add_help = False, usage = '\npython3 seq_select.py -i [input.fasta] -l [list.txt] -o [output.fasta] \npython3 seq_select.py --input [input.fasta] --list [list.txt] --output [output.folder]')
required = parser.add_argument_group('必选项')
optional = parser.add_argument_group('可选项')
required.add_argument('-i', '--input', metavar = '[input.fasta]', help = '输入文件，fasta 格式', required = True)
required.add_argument('-l', '--list', metavar = '[list]', help = '输入文件，txt 格式', required = True)
required.add_argument('-o', '--output', metavar = '[output.folder]', help = '输出文件，fasta 格式', required = True)
#required.add_argument('-l', '--list', metavar = '[list]', help = '记录“序列所在ID/序列起始位置/序列终止位置/以 tab 作为分隔', required = True)
optional.add_argument('-h', '--help', action = 'help', help = '帮助信息')
args = parser.parse_args()

#读取列表文件
list_dict = []
rep5_common = []
rep10_common = []
rep15_common = []
rep20_common = []
rep25_common = []
rep50_common = []
with open(args.list, 'r') as list_file:
    for line in list_file:
        if line.strip():
            line = line.strip().split('\t')
            if int(line[1]) >= 5:
                rep5_common.append(line[0])
            if int(line[1]) >= 10:
                rep10_common.append(line[0])
            if int(line[1]) >= 15:
                rep15_common.append(line[0])
            if int(line[1]) >= 20:
                rep20_common.append(line[0])
            if int(line[1]) >= 25:
                rep25_common.append(line[0])
            if int(line[1]) >= 50:
                rep50_common.append(line[0])
list_file.close()

#**************************************************
#produce repeat seq
rep5_seq = []
rep10_seq = []
rep15_seq = []
rep20_seq = []
rep25_seq = []
rep50_seq = []
for rec in SeqIO.parse(args.input,"fasta"):
    if rec.id in rep5_common:
        rep5_seq.append(rec)
    if rec.id in rep10_common:
        rep10_seq.append(rec)
    if rec.id in rep15_common:
        rep15_seq.append(rec)
    if rec.id in rep20_common:
        rep20_seq.append(rec)
    if rec.id in rep25_common:
        rep25_seq.append(rec)
    if rec.id in rep50_common:
        rep50_seq.append(rec)

result5 = str(args.output+'/'+'rep_5.fa')
result10 = str(args.output+'/'+'rep_10.fa')
result15 = str(args.output+'/'+'rep_15.fa')
result20 = str(args.output+'/'+'rep_20.fa')
result25 = str(args.output+'/'+'rep_25.fa')
result50 = str(args.output+'/'+'rep_50.fa')
SeqIO.write(rep5_seq, result5, "fasta")
SeqIO.write(rep10_seq, result10, "fasta")
SeqIO.write(rep15_seq, result15, "fasta")
SeqIO.write(rep20_seq, result20, "fasta")
SeqIO.write(rep25_seq, result25, "fasta")
SeqIO.write(rep50_seq, result50, "fasta")

