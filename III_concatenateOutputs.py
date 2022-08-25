import argparse
import json
import os
import random
from pprint import pprint

import pandas as pd
from tqdm import trange

from debruijn import read_reads
from IV_sortOutputs import findSupportReadScore

def checkSubSequence(contig, output):
    for o in output:
        if contig in o:
            return False
    return True


def get_args():
    parser = argparse.ArgumentParser()
    # start
    parser.add_argument('-froot', type=str)
    args = parser.parse_args()
    return args



if __name__ == '__main__':
    args = get_args()
    froot = args.froot
    contigs = read_reads(f'{froot}/{froot}.fasta')
    with open(f'{froot}/setting.json') as f:
        setting = json.load(f)
    print(setting)
    filePath = setting['source']
    score_cut = setting['score_cut']
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
                if temp['DENOVO'][i] not in sequences_scores.keys():
                    sequences_scores[temp['DENOVO'][i]] = temp['Score'][i]
                else:
                    sequences_scores[temp['DENOVO'][i]] = temp['Score'][i] + sequences_scores[temp['DENOVO'][i]]

    print('max length before concat: ', len(max(contigs, key=lambda x: len(x))))

    print(len(contigs))
    ks = [i for i in range(20,5,-1)]
    for k in ks:
        nonconcatenatable = [item for item in contigs if len(item) <= k]
        for index in trange(100000):
            concat = False
            concatednatable = [item for item in contigs if item not in nonconcatenatable]
            if len(concatednatable) == 0:
                break
            current = random.choice(concatednatable)
            for i in range(len(contigs)):
                current_head_to_contig_tail_concat = False
                current_tail_to_contig_head_concat = False
                contig = contigs[i]
                if contig == current or len(contig) <= k:
                    continue
                current_head = current[:k]
                contig_tail = contig[len(contig) - k:]
                current_tail = current[len(current)-k:]
                contig_head = contig[:k]
                if current_head == contig_tail: # contig的尾与current的头k个相等
                    contigs[i] = contig + current[k:]
                    concat = True
                    continue
                if current_tail == contig_head: # contig的头与current的尾k个相等
                    contigs[i] = current + contig[k:]
                    concat = True
                    continue
                try:
                    for j in range(k-1):
                        shuffled_current_head = list(current_head)
                        shuffled_current_head[j],shuffled_current_head[j+1] = shuffled_current_head[j+1],shuffled_current_head[j]
                        shuffled_current_head = ''.join(shuffled_current_head)
                        if shuffled_current_head == contig_tail:
                            result_1 = contig + current[k:]
                            result_2 = contig[:len(contig)-k] + current
                            score_1 = findSupportReadScore(result_1,sequences_scores)
                            score_2 = findSupportReadScore(result_2,sequences_scores)
                            if score_1 > score_2:
                                result = result_1
                            else:
                                result = result_2
                            current_head_to_contig_tail_concat = True
                            break
                    if current_head_to_contig_tail_concat:
                        print(current,' ',contig)
                        print(result)
                        contigs[i] = result
                        concat = True
                        continue

                    for j in range(k-1):
                        shuffled_current_tail = list(current_tail)
                        shuffled_current_tail[j],shuffled_current_tail[j+1] = shuffled_current_tail[j+1],shuffled_current_tail[j]
                        shuffled_current_tail = ''.join(shuffled_current_tail)
                        if shuffled_current_tail == contig_head:
                            result_1 = current + contig[k:]
                            result_2 = current[:len(current)-k] + contig
                            score_1 = findSupportReadScore(result_1,sequences_scores)
                            score_2 = findSupportReadScore(result_2,sequences_scores)
                            if score_1 > score_2:
                                result = result_1
                            else:
                                result = result_2
                            current_tail_to_contig_head_concat = True
                            break

                    if current_tail_to_contig_head_concat:
                        print(current,' ',contig)
                        print(result)
                        contigs[i] = result
                        concat = True
                        continue
                except:
                    continue

            if concat:
                contigs.remove(current)
            else:
                nonconcatenatable.append(current)
        print('max length after concat(k = {}): '.format(k), len(max(contigs, key=lambda x: len(x))))
        print(len(contigs))

    outputs = []
    for contig in contigs:
        if contig not in outputs and checkSubSequence(contig, outputs):
            outputs.append(contig)
    print(len(outputs))
    print('max length after concat: ', len(max(outputs, key=lambda x: len(x))))

    k = setting['k_upperlimit']
    outFile = open(f'{froot}/{froot}_concatenated.fasta', mode='a+')
    for i in range(len(outputs)):
        outFile.writelines('>SEQUENCE_{}_{}mer\n{}\n'.format(i, k, outputs[i]))
    outFile.close()
