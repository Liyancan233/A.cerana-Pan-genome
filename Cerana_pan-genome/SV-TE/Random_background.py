#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#_s要统计的字符
#file 文件路径和文件名
#返回统计的数量
#导入模块，初始传递命令、变量等
import argparse
import sys
import re
parser = argparse.ArgumentParser(description = '\n该脚本用于获取随即背景的TE数量', add_help = False, usage = '\n python this scrip.py -l [unique_TEname]  #-o [output.file]  -i [bootstrapN.TMP.txt]\npython3 seq_select.py --input [TE_name_TMP.txt] --output [output.fasta] --list [bootstrap_result]')
required = parser.add_argument_group('必选项')
optional = parser.add_argument_group('可选项')
required.add_argument('-i', '--input', metavar = '[bootstrapN.TMP.txt]', help = '输入文件，TXT 格式', required = True)
#required.add_argument('-o', '--output', metavar = '[output.file]', help = '输出文件，TXT 格式', required = True)
required.add_argument('-l', '--list', metavar = '[unique_TEname.list]', help = '以 tab 作为分隔TXT', required = True)
optional.add_argument('-h', '--help', action = 'help', help = '帮助信息')
args = parser.parse_args()

def cul_str(file,sing):
    with open(args.input, 'r') as f1:
        message=''
        for line in f1:
            message+=line.rstrip()

    counter=len(re.findall(sing,message))
    return counter
#print(len(counter))
#读取列表文件
#print(args.input)

with open(args.list, 'r') as list_file:
    for line in list_file:
        line = line.strip()
        #print(line)
        num = cul_str(args.input,line)
        print(int(num))

