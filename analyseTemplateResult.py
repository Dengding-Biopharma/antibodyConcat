import copy
from pprint import pprint

import pandas as pd
from tqdm import trange


def findOverlap(intervals):
    if len(intervals) <= 1:
        return None
    v1 = max([pair[0] for pair in intervals])
    v2 = min([pair[1] for pair in intervals])
    if v1 > v2:
        return None
    return [v1,v2]


froot = 'avastin_5-10mer_0.6_2'
df = pd.read_csv(f'{froot}/{froot}_blasthomoTemplate.m8', delimiter='\t', header=None)

df = df[df[2] >= 80]

df = df.reset_index(drop=True)

coverage_record = {}
for i in trange(len(df)):
    sequence_label = df[0][i]
    protein_label = df[1][i]
    if sequence_label not in coverage_record.keys():  # new sequence
        protein_record = {}
        coverage_record[sequence_label] = 0
    right = df[9][i]
    left = df[8][i]
    coverage = right - left + 1
    coverage_record[sequence_label] += coverage
    if protein_label not in protein_record.keys():
        protein_record[protein_label] = [[left, right]]
    else:
        protein_record[protein_label].append([left,right])
        overlap = findOverlap(protein_record[protein_label])
        if overlap:
            coverage_record[sequence_label] -= (overlap[1] - overlap[0])
# coverage_record = sorted(coverage_record.items(),key=lambda x:x[1],reverse=True)
print(coverage_record)

outfile = pd.DataFrame()
outfile['sequence ID'] = coverage_record.keys()
outfile['coverage Sum'] = coverage_record.values()
outfile.to_csv(f'{froot}/{froot}_templateCoverageResult.tsv',sep='\t',index=False)

