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
from III_sortOutputs import findSupportReadScore


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
    # parser.add_argument('-rapsearch_path',type=str,required=True)
    # parser.add_argument('-prerapsearch_path', type=str, required=True)
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
                id = line.rstrip()[1:].replace(' ','_')
                contig = lines[i + 1]
                dic[id] = contig.rstrip()
        return dic
    else:
        for i in range(len(lines)):
            line = lines[i]
            if line[0] == '>':
                id = line.rstrip()[1:].replace(' ','_')
                contig = lines[i + 1]
                dic[id] = contig.rstrip()
        return dic


if __name__ == '__main__':
    args = get_args()
    froot = args.froot
    template = args.template
    best_template = open(f'{froot}/best_templates.fasta','w')
    candidates_templates = read_fasta('templates/alpaca.fasta')
    candidates_templates_ann = read_ann('templates/alpaca.ann')
    print(candidates_templates)
    print(candidates_templates_ann)

    template_name = 'templates/alpaca.fasta'
    os.system(f'prerapsearch -d {template_name} -n {froot}/temp-db')
    os.system(f'rapsearch -q {froot}/contigs_sorted.fasta -d {froot}/temp-db -o {froot}/multi_rapsearch_outputs -z 6')
    os.system(
        f'python processRapsearchM8.py -input {froot}/multi_rapsearch_outputs.m8 -output {froot}/multi_rapsearch_outputs_refactor.m8')

    df = pd.read_csv(f'{froot}/multi_rapsearch_outputs_refactor.m8', delimiter='\t', header=None)
    print(df)
    quit()
    for candidates_template in candidates_templates:

        keys = list(region_sequence_dic.keys())
        template_name = f'{froot}/region_temp.fasta'
        with open(template_name,'w') as f:
            for index in trange(len(keys)):
                region_sequence_key = keys[index]
                region_sequence = region_sequence_dic[region_sequence_key]
                f.write(f'>{region_sequence_key}\n{region_sequence}')
        os.system(f'prerapsearch -d {template_name} -n {froot}/temp-db')
        ## contig matching
        os.system(f'rapsearch -q {froot}/contigs_sorted.fasta -d {froot}/temp-db -o {froot}/region_rapsearch_outputs -z 6')
        ## denovo matching
        # os.system(
        #     f'rapsearch -q {froot}/input_reads.fasta -d {froot}/temp-db -o {froot}/region_rapsearch_outputs -z 6')
        os.system(
            f'python processRapsearchM8.py -input {froot}/region_rapsearch_outputs.m8 -output {froot}/region_rapsearch_outputs_refactor.m8')

        df = pd.read_csv(f'{froot}/region_rapsearch_outputs_refactor.m8', delimiter='\t', header=None)
        df = df[df[2] >= 80]
        df = df.reset_index(drop=True)
        for index in trange(len(keys)):
            region_sequence_key = keys[index]
            region_sequence = region_sequence_dic[region_sequence_key]
            try:
                sub_df = df[df[1] == region_sequence_key]
            except:
                continue
            dfList = sub_df.values
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
            # if coverage > 0.97 or index == (len(keys) - 1):
            #     find = True
            region_sequence_coverage_dic = dict(sorted(region_sequence_coverage_dic.items(), key=lambda item: item[1],reverse=True))
            best_template_id = list(region_sequence_coverage_dic.items())[0][0]
            best_coverage = list(region_sequence_coverage_dic.items())[0][1]
            best_template.write(f'>{best_template_id}_{region}_{best_coverage}\n{region_sequence_dic[best_template_id]}\n')