import argparse
import os

import pandas as pd

from generateTemplatesBlastReport import read_fasta


def get_args():
    parser = argparse.ArgumentParser()
    # start
    parser.add_argument('-froot', type=str, required=True)
    # parser.add_argument('-source', type=str, required=True)
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
    best_fragments = read_fasta(f'{froot}/avastin_best_light_fragments.fasta')
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
            if position == 4:
                quit()
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
                step = len(candidate_bases)
                for i in range(num_candidates_letters-1):
                    candidate_bases = candidate_bases + candidate_bases
                for i in range(len(candidate_letters)):
                    candidate_letter = candidate_letters[i]
                    for j in range(i,len(candidate_bases),step):
                        print(i,j)
                        candidate_bases[j].append(candidate_letter)
            print(candidate_bases)
        print(candidate_bases)
        quit()
