import sys
from itertools import islice

#user gaui python this script.py [1]Genename_list.txt [2]go_annot.txt
#[1]Genename_list.txt  formate one clumn; gene list nust start with second line ,frist line is title
#[2] output of R script by interproscan.tsv
 

infile_1 = sys.argv[1]  #gene name
infile_2 = sys.argv[2] #go_annot.txt

f_write = open("./Genename_to_GO.txt", 'w')
with open (infile_1) as f: #Genename  
    for line in islice(f, 1, None):
        lin = line.strip().split('\t')
        with open (infile_2) as g: 
            for lina in g:
                linaa = lina.strip().split('\t')
                gename = linaa[0].split('-mRNA')
                if lin[0] == gename[0]:
                    f_write.write(lin[0] + '\t' + linaa[1] + '\t' + linaa[2] + '\t' + linaa[3] + '\t' + linaa[4] + '\n')
f.close()
g.close()
#f_write.close()
