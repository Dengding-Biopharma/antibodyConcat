import argparse
import json
import os

import pandas as pd
from tqdm import trange

from debruijn import read_reads

def findSupportReadScore(contig,score_table):
    score = 0
    for read in score_table.keys():
        if read in contig:
            score += score_table[read]
    return score

def get_args():
    parser = argparse.ArgumentParser()
    # start
    parser.add_argument('-froot', type=str)
    args = parser.parse_args()
    return args



if __name__ == '__main__':
    args = get_args()
    froot = args.froot
    f = open(f'{froot}/setting.json')
    setting = json.load(f)
    print(setting)
    filePath = setting['source']
    score_cut = setting['score_cut']
    sequences_scores = dict()
    for root, dir, files in os.walk(filePath):
        root = root + '/'
        for file in files:
            filename = root + file
            data = pd.read_csv(filename, delimiter='\t')
            temp = data[data['Score'] >= score_cut]
            temp = temp[-50 < temp['PPM Difference']]
            temp = temp[temp['PPM Difference'] < 50]
            temp.reset_index(inplace=True)
            for i in range(len(temp)):
                if temp['DENOVO'][i] not in sequences_scores.keys():
                    sequences_scores[temp['DENOVO'][i]] = temp['Score'][i]
                else:
                    sequences_scores[temp['DENOVO'][i]] = temp['Score'][i] + sequences_scores[temp['DENOVO'][i]]

    contigs = read_reads(f'{froot}/{froot}.fasta')
    scores = []

    k=setting['k_upperlimit']
    contigs.sort(key=lambda x:findSupportReadScore(x,score_table=sequences_scores),reverse=True)
    outFile = open(f'{froot}/{froot}_sorted.fasta', mode='a+')
    for i in range(len(contigs)):
        outFile.writelines('>SEQUENCE_{}_{}mer_{}\n{}\n'.format(i,k,round(findSupportReadScore(contigs[i],sequences_scores),2),contigs[i]))
    outFile.close()