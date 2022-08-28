import json
from pprint import pprint

import pandas as pd
from tqdm import trange


class aLine:
    def __init__(self, template):
        self.template = list(template)
        self.length = len(self.template)
        self.mask = list('0' * self.length)

    def get_mask(self):
        return ''.join(self.mask)

    def get_template(self):
        return ''.join(self.template)

    def fillMask(self, intervals):
        for interval in intervals:
            for i in range(interval[0] - 1, interval[1]):
                if self.mask[i] == '0':
                    self.mask[i] = '1'

    def get_coverage(self):
        count = 0
        for i in self.mask:
            if i != '0':
                count += 1
        return count / self.length


def read_fasta(path,species=None):
    f = open(path, 'r')
    lines = f.readlines()
    f.close()
    dic = {}
    if species:
        for i in range(len(lines)):
            line = lines[i]
            if line[0] == '>' and species in line:
                id = line.split(' ')[0][1:].rstrip()
                contig = lines[i + 1]
                dic[id] = contig.rstrip()
        return dic
    else:
        for i in range(len(lines)):
            line = lines[i]
            if line[0] == '>':
                id = line.split(' ')[0][1:].rstrip()
                contig = lines[i + 1]
                dic[id] = contig.rstrip()
        return dic

if __name__ == '__main__':

    template_name = 'templates/homo_template.fasta'
    froot = 'avastin_5-10mer_0.6_2'
    contig_filepath = f'{froot}/{froot}_modified_sorted.fasta'
    df = pd.read_csv(f'{froot}/{froot}_blasthomoTemplate.m8', delimiter='\t', header=None)
    template_dic = read_fasta(template_name,'Homo')
    contig_dic = read_fasta(contig_filepath)
    templates = list(template_dic.keys())


    df = df[df[2] >= 80]
    df = df.sort_values(by=1)
    df = df.reset_index(drop=True)

    dfList = df.values
    sequence_template_id_pair_dic={}
    for item in dfList:
        template_id = item[:2][1]
        label = item[:2][0] + '+' + template_id
        value_list = list(item[2:])
        sequence_template_id_pair_dic[label] = value_list


    previous_protein = None
    protein_sequences_records = {}
    for i in trange(len(df)):
        current_protein = df[1][i]
        current_contig_id = df[0][i]
        right = df[9][i]
        left = df[8][i]
        if current_protein not in protein_sequences_records.keys():
            protein_sequences_records[current_protein] = {}
        protein_sequences_records[current_protein][current_contig_id] = [left, right]

    coverage_record = {}
    keys = list(protein_sequences_records.keys())
    for i in trange(len(keys)):
        key = keys[i]
        if key not in templates:
            continue
        protein_sequences_records[key] = dict(sorted(protein_sequences_records[key].items(), key=lambda x: x[1]))
        sequence_interval_records = protein_sequences_records[key]
        intervals = []
        for sub_key in sequence_interval_records.keys():
            intervals.append(sequence_interval_records[sub_key])
        template = template_dic[key]
        line = aLine(template)
        line.fillMask(intervals)
        coverage_record[key] = line.get_coverage()

    coverage_record = dict(sorted(coverage_record.items(), key=lambda x: x[1], reverse=True))

    keys = list(coverage_record.keys())

    with open(f'{froot}/{froot}_templateCoverageResult.json', 'w') as fw:
        for i in trange(len(keys)):
            template_id = keys[i]
            json_block = dict()
            json_block['Template ID'] = template_id
            json_block['Template Sequence'] = template_dic[template_id]
            json_block['Coverage Score'] = coverage_record[template_id]
            matched_contigs = []
            for sub_key in protein_sequences_records[template_id]:
                contig_info = {}
                contig_id = sub_key
                label = contig_id + '+' + template_id
                contig_array = sequence_template_id_pair_dic[label]
                contig_info['Contig ID'] = contig_id
                contig_info['Contig Sequence'] = contig_dic[contig_id]
                contig_info['Identity'] = contig_array[0]
                contig_info['Alignment Length'] = contig_array[1]
                contig_info['Mismatches'] = contig_array[2]
                contig_info['Gap_openings'] = contig_array[3]
                contig_info['Contig Interval'] = str([contig_array[4],contig_array[5]])
                contig_info['Template Interval'] = str([contig_array[6],contig_array[7]])
                contig_info['E Value'] = contig_array[8]
                contig_info['Bit Score'] = contig_array[9]
                matched_contigs += [contig_info]
            json_block['Matched Contigs'] = matched_contigs
            json.dump(json_block, fw, indent=4)
            fw.write('\n')
