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
        df = df[df[2]>90]
        dfList = df.values
        candidate_fragments_dic = {}
        for item in dfList:
            fragment_id = item[1]
            candidate_fragments_dic[fragment_id] = [[item[6]-1,item[7]],[item[8]-1,item[9]]]

        for candidate_fragment in candidate_fragments_dic.keys():
            value = candidate_fragments_dic[candidate_fragment]
            if (value[0][1] - value[0][0]) == (value[1][1] - value[1][0]): # 长度匹配上
                if len(best_fragments[candidate_fragment]) == (value[0][1] - value[0][0]): # 在base上与一块区域重合，包含关系，现在只是continue
                    continue
                elif value[0][1] == len(base): #在tail上match到了，往后拼接
                    print(best_fragments[candidate_fragment])

        quit()