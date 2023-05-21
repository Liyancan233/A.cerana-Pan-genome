import os
import collections
import argparse
from Bio.Seq import Seq
from Bio import SeqIO
#define input_redundant and output_ redundant same name
hit_contigs = []
redundant = []
def best_hit_alignment(alignment_path, identity, coverage):
    f_write = open("./tmp.txt", 'w')
    f = open(alignment_path, 'r')
    P = 1
    #record = 2
    #hit_contigs = set() #use to fill redudant ID
    for line in f.readlines():
        line = line.split('\n')[0]
        line = line.split('\t')
        #line=[(s1)11(e1)2097(s2)2087(e2)1(al-1)2087(al-2)2087(identy)99.95(1-leng)2099(2-leng)2087(1-cov)99.43(2-cov)100.00(1ID)FuJian_LongYan_357357k141_4843(2ID)FuJian_LongYan_358358k141_136 [IDENTITY]]
        cov1 = float(line[9])
        cov2 = float(line[10])
        if line[12] and line[11] not in hit_contigs:
            if float(line[6]) >= float(identity) and (cov2 >= float(coverage)*100 or cov1 >= float(coverage)*100):
                if float(line[7]) >= float(line[8]):
                    hit_contigs.append(line[12])
                    redundant.append(line[11])
                    f_write.write(line[12] + '\t' + line[11] + '\n')
                else:
                    hit_contigs.append(line[11])
                    redundant.append(line[12])
                    f_write.write(line[11] + '\t' + line[12] + '\n')
                    #for num in range(len(line) - 1):
                        #f_write.write(line[num] + '\t')
                    #f_write.write(line[len(line) - 1] + '\n')
    f.close()
    f_write.close()
def extract_ID(input_ref, input_qry, input_oldredundant, output1, output_newredundant):
    newfa2 = []
    newfa3 = []
    ids_ref = []
    ids_qry = []
#******************************************************************************************************
    #produce new_redundant list
    new_redundant = []
    old_redundant = []
    nonchanged_redu = []
    nove_rep = []

    for rec in SeqIO.parse(input_oldredundant,"fasta"):
        old_redundant.append(rec.id)
        if rec.id not in hit_contigs:
            new_redundant.append(rec.id)
        #newfa4.append(rec)
    for red in redundant:
        if red not in new_redundant:
            new_redundant.append(red)
        if red not in old_redundant:
            nove_rep.append(red)
#*****************************************************************************************************
    #produce 1,unredundant fasta 2,unique fasta file
    #ref
    for rec in SeqIO.parse(input_ref,"fasta"):
        if rec.id not in hit_contigs:
            ids_ref.append(rec.id)
    #qry
    for red in SeqIO.parse(input_qry,"fasta"):
        if red.id not in hit_contigs:
            ids_qry.append(red.id)
#*****************************************************************************************************
    #print(ids)_ref
    for rec in SeqIO.parse(input_ref,"fasta"):
        if rec.id in ids_ref:
            newfa2.append(rec)
        if rec.id in new_redundant:
            newfa3.append(rec)
    #print(ids)_qry
    for red in SeqIO.parse(input_qry,"fasta"):
        if red.id in ids_qry:
            newfa2.append(red)
        if red.id in new_redundant:
            newfa3.append(red)
    #print(newfa2,newfa3,newfa4)
    SeqIO.write(newfa2, output1, "fasta")
    SeqIO.write(newfa3, output_newredundant, "fasta")
#***********************************************************************************************************
    #produce new_repeat list
    hit_rep = []
    redu_rep = []
    repeat_write = open("./new_repeat.txt", 'w')
    #f = open("./repeat.txt", 'r')
    #rea = open("./tmp.txt", 'r')
    #f = open("./repeat.txt", 'r+')
    #nove_seq input rep 1
    for rep in nove_rep:
        if rep not in hit_rep:
            repeat_write.write(rep + '\t' + str(1) + '\n')
            hit_rep.append(rep)
    f = open("./repeat.txt", 'r')
    #rea = open("./tmp.txt", 'r')
    for line in f.readlines():
        line = line.split('\n')[0]
        line = line.split('\t')
        #produce nonchanged 
        if (line[0] not in hit_contigs) and (line[0] not in redundant) and (line[0] not in hit_rep):
        #if line[0] not in hit_contigs and line[0] not in redundant and (line[0] not in hit_rep):
            #print(line[0])
            repeat_write.write(line[0] + '\t' + str(line[1]) + '\n')
            hit_rep.append(line[0])
    #produce changing
    f = open("./repeat.txt", 'r')
    for line in f.readlines():
        line = line.split('\n')[0]
        line = line.split('\t')
        rea = open("./tmp.txt", 'r')
        for lina in rea.readlines():
            lina = lina.split('\n')[0]
            lina = lina.split('\t')
            if lina[1] not in redu_rep:
                redu_rep.append(lina[1])
            if lina[0] == line[0] and (lina[0] not in hit_rep) and (lina[1] not in hit_rep):
                repeat_write.write(lina[1] + '\t' + str(int(line[1]) + 1) + '\n')
                hit_rep.append(lina[0])
                hit_rep.append(lina[1])
    #produce changed
    f = open("./repeat.txt", 'r')
    for line in f.readlines():
        line = line.split('\n')[0]
        line = line.split('\t')
        if (line[0] in redu_rep) and (line[0] not in hit_rep):
            repeat_write.write(line[0] + '\t' + str(int(line[1]) + 1) + '\n')
            hit_rep.append(lina[0])
''' #test
    test = open("./test_redu.txt", 'w')
    for lya in new_redundant:
        test.write(lya + '\n')
    allseq = open("./all_redu.txt", 'w') 
    for lyb in ids_ref:
        allseq.write(lyb + '\n')
    for lyd in ids_qry:
        allseq.write(lyd + '\n')
'''
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Filter out alignment results.")
    parser.add_argument("--alignment_path", help="path of coords", required=True, default=None)
    parser.add_argument("--input_ref", help="path of input ref.fa", required=True, default=None)
    parser.add_argument("--input_oldredundant", help="path of input old redundant", required=True, default=None)
    #parser.add_argument("--input_oldrepeat", help="path of input old repeat list ", required=True, default=None)
    parser.add_argument("--input_qry", help="path of input qry.fa", required=True, default=None)
    #parser.add_argument("--REP_bed", help="folder of REP_contigs.bed", required=True, default=None)
    #parser.add_argument("--distance", help="the distance between the placement locations of two contigs",required=False, default=2000)
    parser.add_argument("--identity", help="identity cutoff", required=False, default=0.90)
    parser.add_argument("--coverage", help="coverage cutoff", required=False, default=0.90)
    parser.add_argument("--output1", help="path of nonredundant output", required=True, default=None)
    #parser.add_argument("--output2", help="path of common seq output", required=True, default=None)
    parser.add_argument("--output_newredundant", help="path of new redundant", required=True, default=None)
    FLAGS = parser.parse_args()
    best_hit_alignment(alignment_path=FLAGS.alignment_path, identity=FLAGS.identity, coverage=FLAGS.coverage)
    #contig_pos = placement_pos(BEP_bed=FLAGS.BEP_bed, LEP_bed=FLAGS.LEP_bed, REP_bed=FLAGS.REP_bed)
    extract_ID(input_ref=FLAGS.input_ref, input_oldredundant=FLAGS.input_oldredundant, input_qry=FLAGS.input_qry, output1=FLAGS.output1, output_newredundant=FLAGS.output_newredundant)

