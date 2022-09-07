import argparse
import copy
import json
import os

import pandas as pd

from generateTemplatesBlastReport import read_fasta
from IV_sortOutputs import findSupportReadScore


def get_args():
    parser = argparse.ArgumentParser()
    # start
    parser.add_argument('-froot', type=str, required=True)
    args = parser.parse_args()
    return args


class Aline:
    def __init__(self, sequence):
        self.positions = {}
        for i in range(len(sequence)):
            self.positions[i] = [sequence[i]]


if __name__ == '__main__':
    args = get_args()
    froot = args.froot
    f = open(f'{froot}/setting.json')
    setting = json.load(f)
    print(setting)
    filePath = setting['source']
    score_cut = setting['score_cut']
    best_fragments = read_fasta(f'{froot}/avastin_best_light_fragments.fasta')
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
                if temp['DENOVO'][i].replace('I','L') not in sequences_scores.keys():
                    sequences_scores[temp['DENOVO'][i].replace('I','L')] = temp['Score'][i]
                else:
                    sequences_scores[temp['DENOVO'][i].replace('I','L')] = temp['Score'][i] + sequences_scores[temp['DENOVO'][i].replace('I','L')]
    print(best_fragments)

    base = list(best_fragments.values())[0]
    best_fragments.pop(list(best_fragments.keys())[0])
    print(base)

    while len(best_fragments) != 0:
        with open(f'{froot}/query.fasta', 'w') as f:
            f.write(f'>base\n{base}\n')
        with open(f'{froot}/rest.fasta', 'w') as f:
            for fragment_key in best_fragments.keys():
                f.write(f'>{fragment_key}\n{best_fragments[fragment_key]}\n')

        out = f'{froot}/temp_refactor.m8'
        query = f'{froot}/query.fasta'
        os.system(f'prerapsearch -d {froot}/rest.fasta -n {froot}/rest-db')
        os.system(f'rapsearch -q {query} -d {froot}/rest-db -o {froot}/{froot}_temp')
        os.system(f'python processRapsearchM8.py -input {froot}/{froot}_temp.m8 -output {out}')
        df = pd.read_csv(out, delimiter='\t', header=None)
        df = df[df[2] > 90]
        dfList = df.values
        candidate_fragments_dic = {}
        for item in dfList:
            fragment_id = item[1]
            candidate_fragments_dic[fragment_id] = [[item[6] - 1, item[7]], [item[8] - 1, item[9]]]
        line = Aline(base)
        print(line.positions)
        for candidate_fragment in candidate_fragments_dic.keys():
            value = candidate_fragments_dic[candidate_fragment]
            if (value[0][1] - value[0][0]) == (value[1][1] - value[1][0]):  # 长度匹配上了
                print(value)

                shift = value[0][0]
                for i in range(value[1][0], value[1][1]):
                    if best_fragments[candidate_fragment][i] not in line.positions[i + shift]:
                        line.positions[i + shift].append(best_fragments[candidate_fragment][i])

                if (value[1][1] - value[1][0]) < len(best_fragments[candidate_fragment]):  # fragment可以往后延申
                    for i in range(value[1][1], len(best_fragments[candidate_fragment])):
                        if best_fragments[candidate_fragment][i] not in line.positions[i + shift]:
                            line.positions[i + shift].append(best_fragments[candidate_fragment][i])
        print(line.positions)
        candidate_bases = []
        for position in line.positions.keys():
            candidate_letters = line.positions[position]
            num_candidates_letters = len(candidate_letters)
            if len(candidate_bases) == 0:
                for i in range(num_candidates_letters):
                    candidate_bases.append([candidate_letters[i]])
                continue
            if num_candidates_letters == 1:
                for i in range(len(candidate_bases)):
                    candidate_bases[i].append(candidate_letters[0])
                continue
            if num_candidates_letters > 1:
                candidate_bases_copy = copy.deepcopy(candidate_bases)
                temp = []
                for candidate_letter in candidate_letters:
                    for candidate_base in candidate_bases_copy:
                        candidate_base.append(candidate_letter)
                    temp.extend(candidate_bases_copy)
                    candidate_bases_copy = copy.deepcopy(candidate_bases)
                candidate_bases = temp
        for i in range(len(candidate_bases)):
            candidate_bases[i] = ''.join(candidate_bases[i])


        candidate_bases = sorted(candidate_bases,key=lambda x:findSupportReadScore(x,sequences_scores),reverse=True)
        for candidate_base in candidate_bases:
            print(candidate_base,findSupportReadScore(candidate_bases,sequences_scores))
        base = candidate_bases[0]
        for candidate_fragment in candidate_fragments_dic.keys():
            best_fragments.pop(candidate_fragment)

    print(base)
