import argparse
import copy
import json
import os

import numpy as np
import pandas as pd
from tqdm import trange

from generateTemplatesBlastReport import read_fasta
from III_sortOutputs import findSupportReadScore


def get_args():
    parser = argparse.ArgumentParser()
    # start
    parser.add_argument('-froot', type=str, required=True)
    parser.add_argument('-source', type=str, required=True)
    parser.add_argument('-chain', type=str, required=True)
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
    best_fragments = read_fasta(f'{froot}/{args.source}_best_{args.chain}_fragments.fasta')
    sequences_scores = dict()
    for root, dir, files in os.walk(filePath):
        root = root + '/'
        for file in files:
            filename = root + file
            data = pd.read_csv(filename, delimiter='\t')
            temp = data[data['Score'] >= 0.1]
            temp = temp[-50 < temp['PPM Difference']]
            temp = temp[temp['PPM Difference'] < 50]
            temp.reset_index(inplace=True)
            for i in range(len(temp)):
                if temp['DENOVO'][i].replace('I','L') not in sequences_scores.keys():
                    sequences_scores[temp['DENOVO'][i].replace('I','L')] = temp['Score'][i]
                else:
                    sequences_scores[temp['DENOVO'][i].replace('I','L')] = temp['Score'][i] + sequences_scores[temp['DENOVO'][i].replace('I','L')]

    base = list(best_fragments.values())[0]
    best_fragments.pop(list(best_fragments.keys())[0])
    bases = []
    for index in range(1000):
        current_length = len(best_fragments)
        delete_table = []
        ks = [i for i in range(20, 3, -1)]
        for k in ks:
            for fragment_key in best_fragments.keys():
                fragment = best_fragments[fragment_key]
                if len(fragment) < k:
                    fragment_head = fragment
                else:
                    fragment_head = fragment[:k]
                if len(base) < k:
                    base_tail = base
                else:
                    base_tail = base[len(base)-k:]
                if base_tail == fragment_head:
                    base = base + fragment[k:]
                    delete_table.append(fragment_key)
        for key in delete_table:
            best_fragments.pop(key)
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
        try:
            df = pd.read_csv(out, delimiter='\t', header=None)
            dfList = df.values
            candidate_fragments_dic = {}
            for item in dfList:
                fragment_id = item[1]
                candidate_fragments_dic[fragment_id] = [[item[6] - 1, item[7]], [item[8] - 1, item[9]]]
            line = Aline(base)
            for candidate_fragment in candidate_fragments_dic.keys():
                value = candidate_fragments_dic[candidate_fragment]
                if (value[0][1] - value[0][0]) == (value[1][1] - value[1][0]):  # 长度匹配上了
                    shift = value[0][0]
                    for i in range(value[1][0], value[1][1]):
                        if (i + shift) not in line.positions.keys():
                            line.positions[i + shift] = []
                        if best_fragments[candidate_fragment][i] not in line.positions[i + shift]:
                            line.positions[i + shift].append(best_fragments[candidate_fragment][i])

                    if (value[1][1] - value[1][0]) < len(best_fragments[candidate_fragment]):  # fragment可以往后延申
                        for i in range(value[1][1], len(best_fragments[candidate_fragment])):
                            if (i + shift) not in line.positions.keys():
                                line.positions[i + shift] = []
                            if best_fragments[candidate_fragment][i] not in line.positions[i + shift]:
                                line.positions[i + shift].append(best_fragments[candidate_fragment][i])
            line_keys = list(line.positions.keys())
            position = line_keys[0]
            candidate_letters = line.positions[position]
            num_candidates_letters = len(candidate_letters)
            start_base = []
            for candidate_letter in candidate_letters:
                start_base.append(candidate_letter)
            candidate_bases=start_base
            for position_index in trange(1,len(line_keys)):
                position = line_keys[position_index]
                candidate_letters = line.positions[position]
                num_candidates_letters = len(candidate_letters)
                if num_candidates_letters == 1:
                    for i in range(len(candidate_bases)):
                        candidate_bases[i] = candidate_bases[i] + candidate_letters[0]
                    continue
                if num_candidates_letters > 1:
                    candidate_bases_copy = copy.deepcopy(candidate_bases)
                    # print(num_candidates_letters,len(candidate_bases_copy))
                    candidate_bases = []
                    for candidate_letter in candidate_letters:
                        for candidate_base in candidate_bases_copy:
                            candidate_bases.append(candidate_base+candidate_letter)
                    if len(candidate_bases) > 100000:
                        candidate_bases = sorted(candidate_bases,
                                                 key=lambda x: findSupportReadScore(x, sequences_scores), reverse=True)[:len(candidate_bases)//3]
                    # print('after', len(candidate_bases))
            candidate_bases = sorted(candidate_bases,key=lambda x:findSupportReadScore(x,sequences_scores),reverse=True)
            print('number of candidate bases: ',len(candidate_bases))
            base = candidate_bases[0]
            for candidate_fragment in candidate_fragments_dic.keys():
                best_fragments.pop(candidate_fragment)
        except Exception as e:
            print('current fragments is up to limit!')
            print(base)
            bases.append(base)
            if len(best_fragments) != 0:
                print('left fragments:')
                for best_fragment in best_fragments.keys():
                    print(best_fragments[best_fragment])
                base = list(best_fragments.values())[0]
                best_fragments.pop(list(best_fragments.keys())[0])
                print('new base: ')
                print(base)
            else:
                print('finish!')
                break
    print('##########################final result#########################')
    for base in bases:
        print('fragment: ',base)



