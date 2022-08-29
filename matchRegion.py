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
from generateTemplatesBlastReport import read_fasta
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
    parser.add_argument('-region', type=str,required=True)
    parser.add_argument('-chain', type=str, required=True)
    args = parser.parse_args()
    return args

def NormalizeData(data):
    return (data - np.min(data)) / (np.max(data) - np.min(data))


if __name__ == '__main__':
    args = get_args()
    region = args.region
    chain = args.chain
    froot = args.froot
    template = args.template
    region_file = f'templates/{region}_{template}_{chain}.fasta'
    region_sequence = read_fasta(region_file)
    print(region_sequence)