import sys
from itertools import islice
infile_1 = sys.argv[1] # GFF file
infile_2 = sys.argv[2] # 'choose 'wego' or 'agrigo' 
#user gaid python this.py to extract all gene-name and GO from gff file
#f_write = open("./infile_1.agriGO.txt", 'w')


#******************************************agrigo format*****************************************************************
#******************************************agrigo format*****************************************************************

if infile_2 == 'agrigo':
    with open (infile_1) as f:
        for line in islice(f, 1, None):
            lin = line.strip().split('\t')
            if lin[0] != '###':
                if lin[2] == 'gene':
                    name = lin[8].split(';')[1].split('=')[-1]
                GO = lin[8].split(';')[-2]
                if lin[2] == 'mRNA' and GO.split('=')[0] == 'Ontology_term':
                    #name = lin[8].split(';')[2].split('=')[-1]
                    GOterm = GO.split('=')[1]
                    for num in range(len(GOterm.split(','))):
                        print (name, GOterm.split(',')[num], sep='\t')
                if lin[2] == 'mRNA' and GO.split('=')[0] != 'Ontology_term':
                    print (name, 'None', sep='\t')
    f.close()
    #f_write.close()
    ##output formate
    #augustus_masked-GuiZhou_BiJie_SRR10549553k119_9338-processed-gene-0.0-mRNA-1	GO:0005787
    #augustus_masked-GuiZhou_BiJie_SRR10549553k119_9338-processed-gene-0.0-mRNA-1   GO:0006465
    #augustus_masked-GuiZhou_BiJie_SRR10549553k119_9338-processed-gene-0.0-mRNA-1   GO:0016021
    
#***********************************************wego format***************************************************************
#***********************************************wego format***************************************************************
elif infile_2 == 'wego':

    with open (infile_1) as f:
        for line in islice(f, 1, None):
            lin = line.strip().split('\t')
            if lin[0] != '###':
                if lin[2] == 'gene':
                    name = lin[8].split(';')[1].split('=')[-1]
                GO = lin[8].split(';')[-2]
                if lin[2] == 'mRNA' and GO.split('=')[0] == 'Ontology_term':
                    GOterm = GO.split('=')[1].replace(',', '\t') #将逗号替换为制表符！
                    #f_write.write(name + '\t' + InterPro + '\t' + '\n')
                    print (name, GOterm, sep='\t')
                if lin[2] == 'mRNA' and GO.split('=')[0] != 'Ontology_term':
                    print (name, 'None', sep='\t')
    f.close()
    #f_write.close()
    ##output formate
    #augustus_masked-GuiZhou_BiJie_SRR10549553k119_9338-processed-gene-0.0-mRNA-1	GO:0005787\tGO:0006465\tGO:0016021
else:
    print('incapable of action!! you must choose "wego" or "agrigo" format')