import json
import os
from collections import Counter

import pandas as pd
from tqdm import trange

from debruijn import read_reads


def read_scores(fp):
    f = open(fp, 'r')
    lines = f.readlines()
    f.close()
    scores = []

    for line in lines:
        if line[0] == '>':
            scores.append(float(line.split('_')[-1].strip()))
    return scores


def standard(readInfo):
    output = {}
    output['TITLE'] = readInfo[0]
    output['DENOVO'] = readInfo[1]
    output['Score'] = readInfo[2]
    output['PPM Difference'] = readInfo[3]
    output['Positional Score'] = readInfo[4]
    output['MATCHED'] = readInfo[5]
    return output


df = pd.DataFrame()

froot = 'avastin_5-10mer_0.6_2'
outputfile = f'{froot}/{froot}_modified_sorted.fasta'
settingFile = open(f'{froot}/setting.json', 'r')
setting = json.load(settingFile)
souceFilePath = 'avastin/avastin'
for root, dir, files in os.walk(souceFilePath):
    root = root + '/'
    for file in files:
        filename = root + file
        data = pd.read_csv(filename, delimiter='\t')
        temp = data[data['Score'] >= setting['score_cut']]
        temp = temp[-50 < temp['PPM Difference']]
        temp = temp[temp['PPM Difference'] < 50]
        temp.reset_index(inplace=True, drop=True)
        df = df.append(temp)
df.reset_index(inplace=True, drop=True)

scores = read_scores(outputfile)
contigs = read_reads(outputfile)

data = df.values
with open(f'{froot}/Report.json', 'w') as fw:
    for i in trange(len(contigs)):
        contig = contigs[i]
        json_block = dict()
        support_reads = []
        json_block['index'] = i + 1
        json_block['Contig sequence'] = contig
        json_block['Length'] = len(contig)
        json_block['Score Sum'] = scores[i]
        temp = []
        for j in range(len(df)):
            read = df['DENOVO'][j]
            if read in contig:
                support_reads.append(read)
                readInfo = standard(data[j])
                temp += [readInfo]
        json_block['Supported reads Count(Not Unique)'] = len(support_reads)
        json_block['Supported reads Count(Unique)'] = len(Counter(support_reads))
        json_block['Supported reads Information'] = temp
        json.dump(json_block, fw, indent=4)
        fw.write('\n')
