import os
import collections
import argparse
from Bio.Seq import Seq
from Bio import SeqIO

hit_contigs = []
def best_hit_alignment(alignment_path, identity, coverage, output):
    f_write = open("./blastn_tmp.txt", 'w')
    f = open(alignment_path, 'r')
    #record = 2
    #hit_contigs = set() #use to fill redudant ID
    for line in f.readlines():
        line = line.split('\n')[0]
        line = line.split('\t')

        cov1 = 1.0 * (float(line[5]) - int(line[12])) / float(line[4])
        cov2 = 1.0 * (float(line[5]) - int(line[12])) / float(line[3])
        if line[0] not in hit_contigs:
            if line[0] != line[1]:
                #record = 1
            #else:
                #record = 2
                #hit_contigs.add(line[0])
                #if float(line[2]) <= identity or (float(line[2]) >= identity and (cov1 <= coverage or cov2 <= coverage)):
                if float(line[2]) >= float(identity) and (cov2 >= float(coverage) or cov1 >= float(coverage)):
                    if float(line[3]) >= float(line[4]):
                        hit_contigs.append(line[1])
                    else:
                        hit_contigs.append(line[0])
                #print(hit_contigs)
                    for num in range(len(line) - 1):
                        f_write.write(line[num] + '\t')
                    f_write.write(line[len(line) - 1] + '\n')
    f.close()
def extract_ID(input, hit_contigs, output):
    ids = []
    for rec in SeqIO.parse(input,"fasta"):
        if rec.id not in hit_contigs:
            ids.append(rec.id)
    #print(ids)
    newfa2=[]
    for rec in SeqIO.parse(input,"fasta"):
        if rec.id in ids:
            newfa2.append(rec)
    #print(newfa2)
    SeqIO.write(newfa2, output, "fasta")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Filter out alignment results.")
    parser.add_argument("--alignment_path", help="path of filtered_align.tsv", required=True, default=None)
    parser.add_argument("--input", help="path of input.fa", required=True, default=None)
    #parser.add_argument("--output", help="path of output", required=True, default=None)
    #parser.add_argument("--REP_bed", help="folder of REP_contigs.bed", required=True, default=None)
    #parser.add_argument("--distance", help="the distance between the placement locations of two contigs",required=False, default=2000)
    parser.add_argument("--identity", help="identity cutoff", required=False, default=99)
    parser.add_argument("--coverage", help="coverage cutoff", required=False, default=0.99)
    parser.add_argument("--output", help="path of output.fa", required=True, default=None)
    FLAGS = parser.parse_args()
    best_hit_alignment(alignment_path=FLAGS.alignment_path, identity=FLAGS.identity, coverage=FLAGS.coverage, output=FLAGS.output)
    #contig_pos = placement_pos(BEP_bed=FLAGS.BEP_bed, LEP_bed=FLAGS.LEP_bed, REP_bed=FLAGS.REP_bed)
    extract_ID(input=FLAGS.input, hit_contigs=hit_contigs, output=FLAGS.output)

