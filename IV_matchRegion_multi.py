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
    def __init__(self, template_id, template_sequence, template_region_info):
        self.sequence = template_sequence
        self.id = template_id
        self.match = ['0' for _ in range(len(self.sequence))]
        for region_name, region_interval in template_region_info.items():
            if 'FR' in region_name:
                for i in range(region_interval[0] - 1, region_interval[1] - 2):
                    self.match[i] = 'F'

    def getCoverage(self):
        fr_length = 0
        cover_length = 0
        for position in self.match:
            if position == '1':
                cover_length += 1
                fr_length += 1
            elif position == 'F':
                fr_length += 1

        return cover_length / fr_length


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
    ann = {}

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
            value = [int(value[0]), int(value[1])]
            ann[id][key] = value
    return ann


def get_args():
    parser = argparse.ArgumentParser()
    # start
    parser.add_argument('-froot', type=str, required=True)
    parser.add_argument('-template', type=str, required=True)
    # parser.add_argument('-rapsearch_path',type=str,required=True)
    # parser.add_argument('-prerapsearch_path', type=str, required=True)
    args = parser.parse_args()
    return args


def NormalizeData(data):
    return (data - np.min(data)) / (np.max(data) - np.min(data))


def read_fasta(path, species=None):
    f = open(path, 'r')
    lines = f.readlines()
    f.close()
    dic = {}
    if species:
        for i in range(len(lines)):
            line = lines[i]
            if line[0] == '>' and species in line:
                id = line.rstrip()[1:].replace(' ', '_')
                contig = lines[i + 1]
                dic[id] = contig.rstrip()
        return dic
    else:
        for i in range(len(lines)):
            line = lines[i]
            if line[0] == '>':
                id = line.rstrip()[1:].replace(' ', '_')
                contig = lines[i + 1]
                dic[id] = contig.rstrip()
        return dic


if __name__ == '__main__':
    args = get_args()
    froot = args.froot
    template = args.template
    candidates_templates = read_fasta('templates/alpaca.fasta')
    candidates_templates_ann = read_ann('templates/alpaca.ann')
    # print(candidates_templates)
    # print(candidates_templates_ann)

    template_name = 'templates/alpaca.fasta'
    os.system(f'prerapsearch -d {template_name} -n {froot}/temp-db')
    # denova matching
    os.system(
        f'rapsearch -q {froot}/input_reads.fasta -d {froot}/temp-db -o {froot}/region_rapsearch_outputs -z 6')

    # # contig matching
    # os.system(f'rapsearch -q {froot}/contigs_sorted.fasta -d {froot}/temp-db -o {froot}/multi_rapsearch_outputs -z 6')
    os.system(
        f'python processRapsearchM8.py -input {froot}/multi_rapsearch_outputs.m8 -output {froot}/multi_rapsearch_outputs_refactor.m8')

    df = pd.read_csv(f'{froot}/multi_rapsearch_outputs_refactor.m8', delimiter='\t', header=None)
    df = df[df[2] >= 80]
    df = df.reset_index(drop=True)
    keys = list(candidates_templates_ann.keys())
    Templates = {}
    for key in keys:
        Templates[key] = Template(key, candidates_templates[key], candidates_templates_ann[key])
    region_sequence_coverage_dic = {}
    for i in trange(len(keys)):
        key = keys[i]
        int_key = int(keys[i])
        try:
            sub_df = df[df[1] == int_key]
        except:
            continue
        if sub_df.empty:
            continue
        dfList = sub_df.values
        sequence_template_id_pair_dic = {}
        for item in dfList:
            template_id = str(item[:2][1])
            label = item[:2][0] + '+' + template_id
            value_list = list(item[2:])
            sequence_template_id_pair_dic[label] = value_list

        current_template = Templates[keys[i]]
        for label in sequence_template_id_pair_dic.keys():
            value = sequence_template_id_pair_dic[label]
            if value[5] - value[4] != (value[7] - value[6]):
                continue
            for j in range(value[6] - 1, value[7]):
                if current_template.match[j] == 'F':
                    current_template.match[j] = '1'

        coverage = current_template.getCoverage()
        if coverage >= 0.8:
            region_sequence_coverage_dic[key] = coverage

    region_sequence_coverage_dic = dict(sorted(region_sequence_coverage_dic.items(), key=lambda item: item[1],reverse=True))
    # best_template_id = list(region_sequence_coverage_dic.items())[0][0]
    # best_coverage = list(region_sequence_coverage_dic.items())[0][1]

    keys = list(region_sequence_coverage_dic.keys())
    with open(f'{froot}/best_templates.fasta','w') as f:
        for key in keys:
            f.write(f'>{key}_{region_sequence_coverage_dic[key]}\n{Templates[key].sequence}\n')
