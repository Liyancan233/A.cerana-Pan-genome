import os
import collections
import argparse
from Bio.Seq import Seq
from Bio import SeqIO

hit_contigs = []
def best_hit_alignment(alignment_path, identity, coverage, output):
    f_write = open("./tmp.txt", 'w')
    f = open(alignment_path, 'r')
    #record = 2
    #hit_contigs = set() #use to fill redudant ID
    for line in f.readlines():
        line = line.split('\n')[0]
        line = line.split('\t')

        cov1 = float(line[9])
        cov2 = float(line[10])
        if line[12] and line[11] not in hit_contigs:
            if line[11] != line[12] and float(line[6]) >= float(identity) and (cov2 >= float(coverage)*100 or cov1 >= float(coverage)*100):
                if float(line[7]) >= float(line[8]):
                    hit_contigs.append(line[12])
                else:
                    hit_contigs.append(line[11])
                #print(hit_contigs)
                    for num in range(len(line) - 1):
                        f_write.write(line[num] + '\t')
                    f_write.write(line[len(line) - 1] + '\n')
    f.close()
def extract_ID(input_ref, input_qry, hit_contigs, output):
    #ref
    ids_ref = []
    newfa2=[]
    for rec in SeqIO.parse(input_ref,"fasta"):
        if rec.id not in hit_contigs:
            ids_ref.append(rec.id)
    #print(ids)
    for rec in SeqIO.parse(input_ref,"fasta"):
        if rec.id in ids_ref:
            newfa2.append(rec)
    ids_qry = []
    for red in SeqIO.parse(input_qry,"fasta"):
        if red.id not in hit_contigs:
            ids_qry.append(red.id)
    #print(ids)
    for red in SeqIO.parse(input_qry,"fasta"):
        if red.id in ids_qry:
            newfa2.append(red)
    #print(newfa2)
    SeqIO.write(newfa2, output, "fasta")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Filter out alignment results.")
    parser.add_argument("--alignment_path", help="path of coords", required=True, default=None)
    parser.add_argument("--input_ref", help="path of input ref.fa", required=True, default=None)
    parser.add_argument("--input_qry", help="path of input qry.fa", required=True, default=None)
    parser.add_argument("--identity", help="identity cutoff", required=False, default=90)
    parser.add_argument("--coverage", help="coverage cutoff", required=False, default=0.90)
    parser.add_argument("--output", help="path of output.fa", required=True, default=None)
    FLAGS = parser.parse_args()
    best_hit_alignment(alignment_path=FLAGS.alignment_path, identity=FLAGS.identity, coverage=FLAGS.coverage, output=FLAGS.output)
    extract_ID(input_ref=FLAGS.input_ref, input_qry=FLAGS.input_qry, hit_contigs=hit_contigs, output=FLAGS.output)

