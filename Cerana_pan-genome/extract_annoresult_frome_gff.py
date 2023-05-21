#user gaui python this script.py *.gff(after Maker blastp interproscan anno) > ann.txt 
#two result )(1)-> infile_1.agriGO.txt include gene express and GO (2)-> genename species name
import sys
from itertools import islice
infile_1 = sys.argv[1]
#user gaid python this.py .gff
f_write = open("./infile_1.agriGO.txt", 'w')
with open (infile_1) as f:
    for line in islice(f, 1, None):
        lin = line.strip().split('\t')
        #if lin[0] != '###':
            #name = lin[8].split(';')[1].split('=')[-1]
            #GO = lin[8].split(';')[-2]
            #InterPro = lin[8].split('InterPro:')[-1].split(',')[0]
        if lin[0] != '###' and lin[2] == 'mRNA':
            #Species_name = lin[8].split('(')[-1].split(');')[0]
            InterPro = lin[8].split('Note=')[-1]
            name = lin[8].split(';')[1].split('=')[-1]
            GO = lin[8].split(';')[-2]
            f_write.write(name + '\t' + InterPro + '\t' + GO + '\n')
            leng = InterPro.split(';')
            for num in range(len(leng)):
                retrieval_str = leng[num]
                #print (retrieval_str)
                if retrieval_str == 'Protein of unknown function':
                    print (name,'Unknown',sep='\t')
                    break
                if retrieval_str[len(retrieval_str)-2] == ')':
                #if retrieval_str[-1] == ')':
                    species_name = leng[num].split('(')[-2]
                    print (name,species_name,sep='\t')
                    break
                else:
                    species_name = leng[num].split('(')[-1]
                    print (name,species_name[:-1],sep='\t')
                    break
            #print name, GO,
            #f_write.write('\t' + lin[3] + '\t' + lin[4])
f.close()
f_write.close()