import argparse
import json
import os
from collections import Counter
from pprint import pprint

import numpy as np
import pandas as pd


def get_args():
    parser = argparse.ArgumentParser()
    # start
    parser.add_argument('-source', type=str, required=True)
    parser.add_argument('-score', type=float, required=True)
    parser.add_argument('-t', type=int, required=True)
    parser.add_argument('-kl', type=int, required=True)
    parser.add_argument('-ku', type=int, required=True)
    parser.add_argument('-more', type=int, required=True)
    parser.add_argument('-predfull_path',type=str)
    parser.add_argument('-msSlash_path',type=str)
    parser.add_argument('-spectrum_path',type=str,default=None)
    parser.add_argument('-source_path',type=str,required=True)
    parser.add_argument('-froot',type=str,required=True)
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = get_args()
    more = args.more
    score_cut = args.score
    threshold = args.t
    k_lowerlimit = args.kl
    k_upperlimit = args.ku
    source = f'{args.source}'
    froot = args.froot
    read_path = args.source_path
    setting = {'score_cut': score_cut, 'threshold': threshold, 'k_lowerlimit': k_lowerlimit,
               'k_upperlimit': k_upperlimit,'source':read_path,'more':more,'froot':froot,'source_':source}
    if not more:
        os.mkdir(froot)
        with open(f'{froot}/setting.json', 'w') as fw:
            json.dump(setting, fw, indent=4)
        sequences = []
        f = open(f'{froot}/setting.json')
        setting = json.load(f)
        score_cut = setting['score_cut']
        k_lowerlimit = setting['k_lowerlimit']
        k_upperlimit = setting['k_upperlimit']
        threshold = setting['threshold']
        source = setting['source']
        for root, dir, files in os.walk(source):
            root = root + '/'
            for file in files:
                filename = root + file
                data = pd.read_csv(filename, delimiter='\t')
                temp = data[data['Score'] >= score_cut]
                temp = temp[-50 <= temp['PPM Difference']]
                temp = temp[temp['PPM Difference'] <= 50]
                sequences.extend(temp['DENOVO'].values)

        sequences = list(Counter(sequences).keys())
        with open(f'{froot}/input_reads.fasta', 'w') as fw:
            for i in range(len(sequences)):
                fw.write(f'>input_read{i}\n{sequences[i]}\n')
        quit()
    input_reads = []
    unused_reads =[]
    spectrum_path = args.spectrum_path
    title_denovo_dic = dict()
    for read_filename in os.listdir(read_path):
        read_file = f'{read_path}/{read_filename}'
        data = pd.read_csv(read_file, delimiter='\t')
        for i in range(len(data)):
            title_denovo_dic[data['TITLE'][i]] = [data['DENOVO'][i],data['PPM Difference'][i],data['Score'][i]]

        temp = data[data['Score'] >= score_cut]
        temp = temp[-50 <= temp['PPM Difference']]
        temp = temp[temp['PPM Difference'] <= 50]
        input_reads.extend(temp['DENOVO'].values)

        temp = data[data['Score'] < score_cut]
        temp = temp[temp['Score'] > 0]
        temp = temp[-50 <= temp['PPM Difference']]
        temp = temp[temp['PPM Difference'] <= 50]
        unused_reads.extend(temp['DENOVO'].values)

    unused_reads = list(Counter(unused_reads).keys())
    unused_reads = [x for x in unused_reads if type(x) is str and len(x) > k_lowerlimit]
    print('number of unused reads: ', len(unused_reads))
    df = pd.DataFrame()

    df['Peptide'] = unused_reads
    df['Charge'] = [2 for i in range(len(unused_reads))]
    df['Type'] = ['HCD' for i in range(len(unused_reads))]
    df['NCE'] = [25 for i in range(len(unused_reads))]

    try:
        os.mkdir(froot)
    except:
        pass

    with open(f'{froot}/setting.json', 'w') as fw:
        json.dump(setting, fw, indent=4)
    df.to_csv(f'{froot}/unused_reads.tsv', sep='\t')

    os.system(
        f'python {args.predfull_path} --input {froot}/unused_reads.tsv --model PredFull/pm.h5 --output {froot}/unused_reads_prediction.mgf'
    )
    try:
        os.system(f'touch {froot}/empty.mgf')
    except:
        pass
    df = pd.DataFrame()
    for spectrum_filename in os.listdir(spectrum_path):
        spectrum_file = f'{spectrum_path}/{spectrum_filename}'
        # os.system(
        #     f'./msSLASH/bin/bruteforce  -e {spectrum_file} -l {froot}/unused_reads_prediction.mgf -d {froot}/empty.mgf -o {froot}/msSLASHresult_{spectrum_filename}.tsv'
        # )
        os.system(
            f'./{args.msSlash_path}  -e {spectrum_file} -l {froot}/unused_reads_prediction.mgf -d {froot}/empty.mgf -o {froot}/msSLASHresult_{spectrum_filename}.tsv'
        )
        slashResult = pd.read_csv(f'{froot}/msSLASHresult_{spectrum_filename}.tsv',sep='\t')
        slashResult['DENOVO'] = np.nan
        slashResult['PPM Diff'] = np.nan
        slashResult['Score'] = np.nan
        slashResult.to_csv('test.csv',na_rep=np.nan)
        for i in range(len(slashResult)):
            try:
                data = title_denovo_dic[slashResult['Title'][i]]
                slashResult['DENOVO'][i] = data[0]
                slashResult['PPM Diff'][i] = data[1]
                slashResult['Score'][i] = data[2]
            except:
                continue
        df = df.append(slashResult)
        df.reset_index(drop=True, inplace=True)

    temp = df[df['TopScore'] >= 0.5]
    temp = temp[temp['Score'] >= 0.6]
    temp = temp[temp['PPM Diff'] >= -50]
    temp = temp[temp['PPM Diff'] <= 50]
    unused_reads = temp['TopPep'].values
    input_reads.extend(unused_reads)
    input_reads = list(Counter(input_reads).keys())
    input_reads = [x for x in input_reads if type(x) is str and len(x) > k_lowerlimit]
    with open(f'{froot}/input_reads.fasta', 'w') as fw:
        for i in range(len(input_reads)):
            fw.write(f'>input_read{i}\n{input_reads[i]}\n')

