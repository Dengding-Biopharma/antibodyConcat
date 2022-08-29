import argparse
import json
import os
import re
import uuid
from collections import Counter
from pprint import pprint
import ast
import numpy as np
import subprocess
import pandas as pd
from tqdm import trange
from IV_sortOutputs import findSupportReadScore


class Template:
    def __init__(self, template_id, template_sequence):
        self.sequence = template_sequence
        self.id = template_id
        self.contigArrays = []
        self.different_position = []
        self.letters_errorRate = {}
        for i in range(len(self.sequence)):
            self.letters_errorRate[i] = {}
        self.unusedReads_match = {}
        for i in range(len(self.sequence)):
            self.unusedReads_match[i] = []


class Contig:
    def __init__(self, contig_id, contig_sequence, template_interval, contig_interval):
        self.sequence = contig_sequence
        self.id = contig_id
        self.template_interval = template_interval
        self.contig_interval = contig_interval
        self.rates = {}
        for i in range(len(self.sequence)):
            self.rates[i] = 0


class fillingTemplate:
    def __init__(self, template_sequence):
        self.template_sequence = template_sequence
        self.fill = list(' ' * len(self.template_sequence))

    def fill_match(self, contig):
        contig_sequence = list(contig.sequence[contig.contig_interval[0] - 1:contig.contig_interval[1]])
        self.fill[contig.template_interval[0] - 1:contig.template_interval[1]] = contig_sequence

    def get_match_result(self):
        return ''.join(self.fill)


def checkOverlap(contig_array, contig):
    intervals = []
    for item in contig_array:
        intervals.append(item.template_interval)
    intervals.append(contig.template_interval)
    intervals.sort(key=lambda x: x[0])

    for i in range(len(intervals) - 1):
        if intervals[i][1] > intervals[i + 1][0]:
            return True
    return False



def read_ann(file_path):
    f = open(file_path, 'r')
    lines = f.readlines()
    f.close()
    ann={}

    for i in range(len(lines)):
        line = lines[i]
        if '>' in line:
            id = line.rstrip()[1:]
            ann[id] = {}
            continue
        else:
            fragment = line.split('=')
            key = fragment[0]
            temp = fragment[1].rstrip()
            value = temp.split('-')
            value = [int(value[0])-1,int(value[1])-2]
            ann[id][key] = value
    return ann

def get_args():
    parser = argparse.ArgumentParser()
    # start
    parser.add_argument('-froot', type=str,required=True)
    parser.add_argument('-template', type=str,required=True)
    args = parser.parse_args()
    return args

def NormalizeData(data):
    return (data - np.min(data)) / (np.max(data) - np.min(data))

def read_fasta(path,species=None):
    f = open(path, 'r')
    lines = f.readlines()
    f.close()
    dic = {}
    if species:
        for i in range(len(lines)):
            line = lines[i]
            if line[0] == '>' and species in line:
                id = line.rstrip()[1:]
                contig = lines[i + 1]
                dic[id] = contig.rstrip()
        return dic
    else:
        for i in range(len(lines)):
            line = lines[i]
            if line[0] == '>':
                id = line.rstrip()[1:]
                contig = lines[i + 1]
                dic[id] = contig.rstrip()
        return dic


if __name__ == '__main__':
    args = get_args()
    regions = ['NonConstant','constant',]
    chains = ['Heavy', 'Light']
    froot = args.froot
    template = args.template
    best_template = open(f'{froot}/best_templates','w')
    for region in regions:
        for chain in chains:
            region_file = f'templates/{region}_{template}_{chain}.fasta'
            region_sequence_dic = read_fasta(region_file)
            print(list(region_sequence_dic.items())[0])
            print(list(region_sequence_dic.items())[0][0])
            print(list(region_sequence_dic.items())[0][1])
            quit()
            region_sequence_coverage_dic = {}
            keys = list(region_sequence_dic.keys())
            for index in trange(len(keys)):
                region_sequence_key = keys[index]
                region_sequence = region_sequence_dic[region_sequence_key]
                template_name = f'{froot}/region_temp.fasta'
                with open(template_name,'w') as f:
                    f.write(f'>{region_sequence_key}\n{region_sequence}')
                os.system(f'prerapsearch -d {template_name} -n {froot}/temp-db')
                os.system(f'rapsearch -q {froot}/{froot}_sorted.fasta -d {froot}/temp-db -o {froot}/region_rapsearch_outputs -z 6')
                os.system(
                    f'python processRapsearchM8.py -input {froot}/region_rapsearch_outputs.m8 -output {froot}/region_rapsearch_outputs_refactor.m8')
                try:
                    df = pd.read_csv(f'{froot}/region_rapsearch_outputs_refactor.m8', delimiter='\t', header=None)
                except:
                    continue
                df = df[df[2] >= 80]
                df = df.reset_index(drop=True)
                dfList = df.values
                sequence_template_id_pair_dic = {}
                for item in dfList:
                    template_id = item[:2][1]
                    label = item[:2][0] + '+' + template_id
                    value_list = list(item[2:])
                    sequence_template_id_pair_dic[label] = value_list
                blank_sequence = list('0'*len(region_sequence))
                for label in sequence_template_id_pair_dic.keys():
                    value = sequence_template_id_pair_dic[label]
                    if value[5]-value[4] != (value[7]-value[6]):
                        continue
                    for i in range(value[6]-1,value[7]):
                        blank_sequence[i] = '1'
                coverage = blank_sequence.count('1')/len(blank_sequence)
                region_sequence_coverage_dic[region_sequence_key] = coverage
            region_sequence_coverage_dic = dict(sorted(region_sequence_coverage_dic.items(), key=lambda item: item[1],reverse=True))
            best_template.write(f'>{list(region_sequence_coverage_dic.items())[0][0]}_{region}\n{list(region_sequence_coverage_dic.items())[0][1]}')