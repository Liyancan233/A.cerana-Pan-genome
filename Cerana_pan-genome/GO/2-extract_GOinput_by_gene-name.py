import sys
from itertools import islice

#user gaui python this script.py [1]Genename_list.txt [2]GO background
#[1]Genename_list.txt  formate one clumn; gene list nust frist line is title
#[2] output as for wego or agrigo input
#[1]Genename_list.txt
#Genename
#gene-name-01
infile_1 = sys.argv[1] #gene name
infile_2 = sys.argv[2] #wego/agrigo format background
unique = []
#f_write = open("./Genename_to_GO.txt", 'w')
with open (infile_1) as f: #Genename  
    for line in islice(f, 1, None):
        lin = line.strip().split('\t')
        with open (infile_2) as g: #background
            for lina in g:
                linaa = lina.strip().split('\t')
                if lin[0] == linaa[0] and lina not in unique:
                    print(lina.strip())
                    unique.append(lina)
f.close()
g.close()
#f_write.close()
