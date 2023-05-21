#user gaui python this script.py infile1 infile2 > result.txt 
#两个vcf去除表头的前七列，主要需要1染色体，2位置，该脚本主要是用在SURVING合并多个vcf之后，如果差500bp的情况下可以从delly中把附近的位点取出
import sys
infile_1 = sys.argv[1] ##delly sv result 
infile_2 = sys.argv[2] ##Surving result delly smoove manta combine result
# cat SURVIVOR-Final_merged.vcf | grep -v '^#' | awk 'BEGIN{OFS="\t"; FS="\t"} {print $1,$2,$3,$4,$5,$6,$7}' > 1-7line.txt
#f_write = open("./infile_1.agriGO.txt", 'w')
with open (infile_1) as f:
    for line in f:
        lined = line.strip().split('\t')
        #print(lined[0],lined[1])
        with open (infile_2) as g:
            for lin in g:
                lind = lin.strip().split('\t')
                #print(lind[0],lind[1])
                if lined[0] == lind[0]:
                    if int(lind[1]) in range(int(lined[1])-500,int(lined[1])+500):
                    #if  float(lined[1])-500<=lind[1]<=float(lined[1])+500:
                        #print(lined[0])
                        print(lind[0],lind[1],sep='\t')
f.close()
