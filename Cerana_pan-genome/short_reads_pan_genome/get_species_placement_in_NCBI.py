import sys
from ete3 import NCBITaxa
input_file = sys.argv[1]
output_file = sys.argv[2]
ncbi = NCBITaxa()
f = open("out.txt", "w") 
fw = open(output_file,"w")
with open(input_file,"r") as fr:
    for line in fr:
        if len(line.strip().split('\t')) == 2:
            species_name = line.strip().split('\t')[1]
            gene_name = line.strip().split('\t')[0]
        else:
            species_name = line.strip()
            gene_name = species_name
        name2taxid = ncbi.get_name_translator([species_name])
        for a,b in name2taxid.items():
            lineage = ncbi.get_lineage(b[0])
            names = ncbi.get_taxid_translator(lineage)
            if names is None:
                fw.write(gene_name+"\t"+'none'+'\n')
            else:
                fw.write(gene_name+"\t")
            i = 1
            #fw.write(species_name+"\t")
            for taxid in lineage:   
                if i < len(lineage):
                    fw.write(names[taxid]+",")
                    i = i + 1
                else:
                    fw.write(names[taxid]+"\n")
        print>>f,(species_name + "OK")


f.close()        
fw.close()

