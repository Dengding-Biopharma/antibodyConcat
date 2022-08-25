import argparse
import copy
import json
import os

from collections import Counter
import sys
from pprint import pprint

import pandas as pd
import debruijn as db
from debruijn import read_reads

def getScore(edge_count_table, contig, k):
    score = 0
    for i in range(len(contig) - k):
        score += edge_count_table[contig[i:i + k + 1]]
    return score

def get_args():
    parser = argparse.ArgumentParser()
    # start
    parser.add_argument('-froot', type=str)
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    sequences = []
    args = get_args()
    froot = args.froot
    # froot = 'avastin_5-8mer_0.8_2'
    f = open(f'{froot}/setting.json')
    setting = json.load(f)
    score_cut = setting['score_cut']
    k_lowerlimit = setting['k_lowerlimit']
    k_upperlimit = setting['k_upperlimit']
    threshold = setting['threshold']
    source = setting['source']
    for root, dir, files in os.walk(source):
        root = root + '/'
        for file in files:
            filename = root + file
            data = pd.read_csv(filename, delimiter='\t')
            temp = data[data['Score'] >= score_cut]
            temp = temp[-50<=temp['PPM Difference']]
            temp = temp[temp['PPM Difference']<=50]
            sequences.extend(temp['DENOVO'].values)




    sequences = list(Counter(sequences).keys())
    # sequences = read_reads(f'{froot}/input_reads.fasta')
    print(len(sequences))
    # sequences = ['EVQLVE','QLVAPG','LVESGGAL','LVESGGGL']
    for k in range(k_lowerlimit, k_upperlimit + 1):
        if k <= k_upperlimit - 1:
            g, pull_out_read, branch_kmer, already_pull_out, edge_count_table = db.construct_graph(sequences, k,
                                                                                                   threshold=threshold)
        else:
            g, pull_out_read, branch_kmer, already_pull_out, edge_count_table = db.construct_graph(sequences, k,
                                                                                                   threshold=threshold, final=True)
        sequences = db.output_contigs(g, branch_kmer, already_pull_out)
        sequences.sort(key=lambda x: getScore(edge_count_table, x, k), reverse=True)
        if k == k_upperlimit:
            outFile = open(f'{froot}/{froot}.fasta', mode='a+')
            for i in range(len(sequences)):
                outFile.writelines('>SEQUENCE_{}_{}mer\n{}\n'.format(i, k, sequences[i]))
            outFile.close()
            break
        print('max length: ', len(max(sequences, key=lambda x: len(x))))
        print('number of output for k={}: '.format(k), len(sequences))
        if k <= k_upperlimit - 1:
            sequences.extend(pull_out_read)
            print('number of pull out read: ', len(pull_out_read))

