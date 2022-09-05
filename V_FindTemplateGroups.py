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
import naive_debruijn_graph as naive_db


class Template:
    def __init__(self, template_id, template_sequence, template_type):
        self.sequence = template_sequence
        self.type = template_type
        self.ignore = True if type == 'nc' else False
        self.id = template_id
        self.contigArrays = []
        self.different_position = []
        self.letters_errorRate = {}
        for i in range(len(self.sequence)):
            self.letters_errorRate[i] = {}
        self.unusedReads_match = {}
        for i in range(len(self.sequence)):
            self.unusedReads_match[i] = []
        self.best_fragments = []


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
            value = [int(value[0]) - 1, int(value[1]) - 1]
            ann[id][key] = value
    return ann


def get_args():
    parser = argparse.ArgumentParser()
    # start
    parser.add_argument('-froot', type=str, required=True)
    parser.add_argument('-source', type=str, required=True)
    args = parser.parse_args()
    return args


def NormalizeData(data):
    return (data - np.min(data)) / (np.max(data) - np.min(data))


if __name__ == '__main__':
    args = get_args()
    froot = args.froot
    template_name = f'{froot}/best_templates.fasta'
    annotation_name = f'templates/mAB_database.ann'
    contig_filepath = f'{froot}/{froot}_sorted.fasta'
    settingFile = open(f'{froot}/setting.json', 'r')
    setting = json.load(settingFile)
    sourceFilePath = f'{args.source}/{args.source}'
    DF = pd.DataFrame()
    unused_reads = pd.DataFrame()
    sequences_scores = dict()
    for root, dir, files in os.walk(sourceFilePath):
        root = root + '/'
        for file in files:
            filename = root + file
            data = pd.read_csv(filename, delimiter='\t')

            unused_reads_temp = data[data['Score'] >= 0.5]
            unused_reads_temp = unused_reads_temp[-50 < unused_reads_temp['PPM Difference']]
            unused_reads_temp = unused_reads_temp[unused_reads_temp['PPM Difference'] < 50]
            unused_reads_temp.reset_index(inplace=True, drop=True)
            unused_reads = unused_reads.append(unused_reads_temp)

            temp = data[data['Score'] >= setting['score_cut']]
            temp = temp[-50 < temp['PPM Difference']]
            temp = temp[temp['PPM Difference'] < 50]
            temp.reset_index(inplace=True, drop=True)
            for i in range(len(temp)):
                if temp['DENOVO'][i].replace('I', 'L') not in sequences_scores.keys():
                    sequences_scores[temp['DENOVO'][i].replace('I', 'L')] = temp['Score'][i]
                else:
                    sequences_scores[temp['DENOVO'][i].replace('I', 'L')] = temp['Score'][i] + sequences_scores[
                        temp['DENOVO'][i].replace('I', 'L')]
            DF = DF.append(temp)

    DF.reset_index(inplace=True, drop=True)
    unused_reads.reset_index(inplace=True, drop=True)

    os.system(f'prerapsearch -d {template_name} -n templates/temp-db')
    os.system(f'rapsearch -q {froot}/{froot}_sorted.fasta -d templates/temp-db -o {froot}/rapsearch_outputs -z 6')
    os.system(
        f'python processRapsearchM8.py -input {froot}/rapsearch_outputs.m8 -output {froot}/rapsearch_outputs_refactor.m8')
    df = pd.read_csv(f'{froot}/rapsearch_outputs_refactor.m8', delimiter='\t', header=None)
    # df = pd.read_csv(f'{froot}/{froot}_blasthomoTemplate.m8',delimiter='\t',header=None)
    df = df[df[2] >= 80]
    df = df.sort_values(by=0)
    df = df.reset_index(drop=True)

    template_dic = read_fasta(template_name)
    templates = list(template_dic.keys())
    contig_dic = read_fasta(contig_filepath)
    contigs = list(contig_dic.keys())
    annotation = read_ann(annotation_name)

    dfList = df.values
    sequence_template_id_pair_dic = {}
    for item in dfList:
        template_id = item[:2][1]
        label = item[:2][0] + '+' + template_id
        value_list = list(item[2:])
        sequence_template_id_pair_dic[label] = value_list

    template_contig_group = {}
    while len(contigs) != 0:
        current_contig = contigs[0]
        best_identity = 0
        best_template = None
        longest_length = 0
        for template_id in templates:
            label = current_contig + '+' + template_id
            try:
                identity = sequence_template_id_pair_dic[label][0]
            except:
                continue
            # if identity > 90:
            #     template_length_record[template_id] = len(template_dic[template_id])
            if identity >= best_identity and identity > 90 and len(template_dic[template_id]) > longest_length:
                longest_length = len(template_dic[template_id])
                best_identity = identity
                best_template = template_id
        # if len(template_length_record) == 0:
        if not best_template:
            contigs.remove(current_contig)
            continue
        # candidate_templates = sorted(list(template_length_record.items()), key=lambda x: x[1], reverse=True)
        # best_template = candidate_templates[0][0]
        template_contig_group[best_template] = [current_contig]
        contigs.remove(current_contig)

        remove = []
        for contig in contigs:
            label = contig + '+' + best_template
            try:
                value = sequence_template_id_pair_dic[label]
            except:
                continue
            if value:
                template_contig_group[best_template] += [contig]
                remove.append(contig)
        for item in remove:
            contigs.remove(item)

    # pprint(template_contig_group)
    # quit()
    report_path = f'{froot}/{froot}_TemplateMatchReport.txt'
    outFile = open(report_path, 'w')
    message = ''
    html_path = f'{froot}/{froot}_TemplateMatchReport.html'
    htmlFile = open(html_path, 'w')
    html = '''<!DOCTYPE html>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/Chart.js/2.9.4/Chart.js"></script>
    <body>
    '''
    reads = DF['DENOVO'].values
    unused_reads = unused_reads['DENOVO'].values
    reads_count = Counter(reads)
    position_scores = DF['Positional Score'].values

    unused_reads = list(Counter(unused_reads))

    with open(f'{froot}/unusedReads.fasta', 'w') as f:
        for i in range(len(unused_reads)):
            f.write('>unused_reads_{}\n'.format(i))
            f.write('{}\n'.format(unused_reads[i]))

    Templates = []
    for template_id in template_contig_group.keys():
        type = 'nc' if 'NonConstant' in template_id else 'c'
        template = Template(template_id, template_dic[template_id].replace('I', 'L'), type)
        Templates.append(template)
        for contig_id in template_contig_group[template_id]:
            label = contig_id + '+' + template_id
            value = sequence_template_id_pair_dic[label]
            contig = Contig(contig_id, contig_dic[contig_id], [value[6], value[7]], [value[4], value[5]])
            if len(template.contigArrays) > 0:
                overlap = True
                for contig_array in template.contigArrays:
                    if not checkOverlap(contig_array, contig):
                        contig_array.append(contig)
                        overlap = False
                        break
                if overlap:
                    template.contigArrays.append([contig])
            else:
                template.contigArrays.append([contig])

        for array_index in range(len(template.contigArrays)):
            contig_array = template.contigArrays[array_index]
            contig_array = sorted(contig_array, key=lambda x: x.template_interval[0])
            template.contigArrays[array_index] = contig_array
            for contig in contig_array:
                for i in range(len(reads)):
                    read = reads[i]
                    if read in contig.sequence:
                        match = re.search(read, contig.sequence)
                        read_positional_scores = ast.literal_eval(position_scores[i])
                        for j in range(match.start(), match.end()):
                            positional_score = 1 - read_positional_scores[j - match.start()]
                            if positional_score == 0:
                                contig.rates[j] += -20
                            else:
                                contig.rates[j] += np.log((1 - read_positional_scores[j - match.start()]))

        minimum_contigs_array = []
        start_contig = template.contigArrays[0][0]
        for i in range(1, len(template.contigArrays)):
            if template.contigArrays[i][0].template_interval[0] < start_contig.template_interval[0] \
                    and \
                    template.contigArrays[i][0].template_interval[1] > start_contig.template_interval[1]:
                start_contig = template.contigArrays[i][0]
        minimum_contigs_array.append(start_contig)
        current_contig = start_contig
        candidate_contigs = []
        while True:
            for array_index in range(len(template.contigArrays)):
                contig_array = template.contigArrays[array_index]
                for contig_index in range(len(contig_array)):
                    contig = contig_array[contig_index]
                    if current_contig.template_interval[0] < contig.template_interval[0] < \
                            current_contig.template_interval[1] < \
                            contig.template_interval[1]:
                        candidate_contigs.append(contig)
            if len(candidate_contigs) == 0:
                for array_index in range(len(template.contigArrays)):
                    contig_array = template.contigArrays[array_index]
                    for contig_index in range(len(contig_array)):
                        contig = contig_array[contig_index]
                        if current_contig.template_interval[1] <= contig.template_interval[0]:
                            candidate_contigs.append(contig)
                if len(candidate_contigs) == 0:
                    break
                candidate_contigs = sorted(candidate_contigs, key=lambda x: x.template_interval[0])
                current_contig = candidate_contigs[0]
                minimum_contigs_array.append(current_contig)
                candidate_contigs = []
            else:
                candidate_contigs = sorted(candidate_contigs, key=lambda x: x.template_interval[1], reverse=True)
                current_contig = candidate_contigs[0]
                minimum_contigs_array.append(current_contig)
                candidate_contigs = []

        # print()
        # print('*' * 500)
        # print(template.sequence)
        message += '*' * 500
        message += '\n'
        message += 'Template ID: {}'.format(template.id)
        message += '\n'
        message += template.sequence
        message += '\n'
        for contig_array in template.contigArrays:
            contig_array = sorted(contig_array, key=lambda x: x.template_interval[0])
            for contig in contig_array:
                if (contig.contig_interval[1] - contig.contig_interval[0]) != (
                        contig.template_interval[1] - contig.template_interval[0]):
                    continue
                template_points = [x for x in range(contig.template_interval[0] - 1, contig.template_interval[1])]
                contig_points = [x for x in range(contig.contig_interval[0] - 1, contig.contig_interval[1])]
                for i in range(len(template_points)):
                    template_point = template_points[i]
                    contig_point = contig_points[i]
                    current_template_position = template.letters_errorRate[template_point]
                    current_contig_letter = contig.sequence[contig_point]
                    if current_contig_letter not in current_template_position.keys():
                        current_template_position[current_contig_letter] = contig.rates[contig_point]
                    else:
                        current_template_position[current_contig_letter] = current_template_position[
                                                                               current_contig_letter] + contig.rates[
                                                                               contig_point]

        position_keys = list(template.letters_errorRate.keys())
        result_sequences = []
        for key in position_keys:
            candidate_letters = template.letters_errorRate[key]
            candidate_letters = dict(sorted(candidate_letters.items(), key=lambda item: item[1]))
            if candidate_letters != {}:
                while len(result_sequences) < len(candidate_letters):
                    result_sequences.append(list(' ' * len(template.sequence)))
                letters = list(candidate_letters.keys())
                for i in range(len(letters)):
                    result_sequences[i][key] = letters[i]
        for sequence in result_sequences:
            # print(''.join(sequence))
            message += ''.join(sequence)
            message += '\n'

        message += '-' * 100 + '\n'
        # print('-' * 100)
        message += 'Unused reads blast result: \n'
        # print('Unused reads blast result: ')
        message += template.sequence + '\n'
        # print(template.sequence)
        with open(f'{froot}/temp.fasta', 'w') as f:
            f.write('>{}\n'.format(template.id))
            f.write(template.sequence)

        out = f'{froot}/{froot}_unusedReadsBlastTemplate_refactor.m8'
        query = f'{froot}/unusedReads.fasta'
        os.system(f'prerapsearch -d {froot}/temp.fasta -n {froot}/temp')
        os.system(f'rapsearch -q {query} -d {froot}/temp -o {froot}/{froot}_unusedReadsBlastTemplate')
        os.system(f'python processRapsearchM8.py -input {froot}/{froot}_unusedReadsBlastTemplate.m8 -output {out}')

        unusedReadsTemplateResults = pd.read_csv(out, delimiter='\t', header=None)
        unusedReadsTemplateResults = unusedReadsTemplateResults[unusedReadsTemplateResults[2] >= 90]
        unusedReadsTemplateResults.reset_index(drop=True, inplace=True)

        unusedReads_dic = read_fasta(query)

        unusedReads_value_dic = {}
        for value in unusedReadsTemplateResults.values:
            unusedReads_value_dic[value[0]] = list(value[1:])
        unused_reads_intervals = {}
        for unusedRead in unusedReads_value_dic.keys():
            item = unusedReads_value_dic[unusedRead]
            read_left = item[5]
            read_right = item[6]
            template_left = item[7]
            template_right = item[8]
            if (read_right - read_left) != (template_right - template_left):
                continue
            unused_reads_intervals[unusedReads_dic[unusedRead]] = [[template_left, template_right],
                                                                   [read_left, read_right]]
            matchedReadSeq = unusedReads_dic[unusedRead][read_left - 1:read_right]
            for i in range(template_left - 1, template_right):
                current_read_letter = matchedReadSeq[i - (template_left - 1)]
                if current_read_letter not in template.unusedReads_match[i]:
                    template.unusedReads_match[i] += [current_read_letter]
        unusedReadsResultSequence = []
        for key in template.unusedReads_match.keys():
            while len(template.unusedReads_match[key]) > len(unusedReadsResultSequence):
                unusedReadsResultSequence.append(list(' ' * len(template.sequence)))
            letters = template.unusedReads_match[key]
            for i in range(len(letters)):
                unusedReadsResultSequence[i][key] = letters[i]
        for sequence in unusedReadsResultSequence:
            # print(''.join(sequence))
            message += ''.join(sequence) + '\n'

        matched_length = 0
        for i in range(len(template.sequence)):
            if template.letters_errorRate[i] or template.unusedReads_match[i]:
                matched_length += 1

        coverage = matched_length / len(template.sequence)
        if coverage < 0.5:
            continue

        merged_result = []
        while len(result_sequences) < max(len(result_sequences), len(unusedReadsResultSequence)):
            result_sequences.append(list(' ' * len(template.sequence)))
        for unusedReadResult in unusedReadsResultSequence:
            for i in range(len(template.sequence)):
                letter = unusedReadResult[i]
                for contig_result in result_sequences:
                    if contig_result[i] == letter:
                        break
                    if contig_result[i] == ' ':
                        contig_result[i] = '<font color="green">{}</font>'.format(letter)
                        break

        for sequence in result_sequences:
            for i in range(len(sequence)):
                if sequence[i] != ' ' and len(sequence[i]) == 1:
                    sequence[i] = '<font color="blue">{}</font>'.format(sequence[i])

        merged_result = result_sequences

        best_result = merged_result[0]
        best_result_fragments = []
        counting = False
        for i in range(len(best_result)):
            best_result_position = best_result[i]
            if not counting and best_result_position != ' ':
                fragment = ''
                counting = True
                fragment += ''.join(c for c in best_result_position if c.isupper())
            elif counting and best_result_position != ' ' and i != (len(best_result) - 1):
                fragment += ''.join(c for c in best_result_position if c.isupper())
            elif (counting and best_result_position == ' ') and i != (len(best_result) - 1):
                best_result_fragments.append(fragment)
                counting = False
            elif i == (len(best_result) - 1):
                fragment += ''.join(c for c in best_result_position if c.isupper())
                best_result_fragments.append(fragment)
                counting = False

        k = 25
        best_contigs = []
        for fragment in best_result_fragments:
            if len(fragment) <= k:
                head = fragment
                tail = fragment
            else:
                head = fragment[:k]
                tail = fragment[len(fragment) - k:]
            with open(f'{froot}/head.fasta', 'w') as f:
                f.write(f'>head\n{head}')
            with open(f'{froot}/tail.fasta', 'w') as f:
                f.write(f'>tail\n{tail}')
            head_out = f'{froot}/{froot}_head_best_contigs_refactor.m8'
            tail_out = f'{froot}/{froot}_tail_best_contigs_refactor.m8'
            head_query = f'{froot}/head.fasta'
            tail_query = f'{froot}/tail.fasta'
            os.system(f'prerapsearch -d {contig_filepath} -n {froot}/contigs')
            os.system(f'rapsearch -q {head_query} -d {froot}/contigs -o {froot}/{froot}_head_best_contigs')
            os.system(f'rapsearch -q {tail_query} -d {froot}/contigs -o {froot}/{froot}_tail_best_contigs')
            os.system(f'python processRapsearchM8.py -input {froot}/{froot}_head_best_contigs.m8 -output {head_out}')
            os.system(f'python processRapsearchM8.py -input {froot}/{froot}_tail_best_contigs.m8 -output {tail_out}')
            try:
                # if not template.ignore:
                template.ignore = False
                head_df = pd.read_csv(head_out, delimiter='\t', header=None)
                head_df = head_df[head_df[2] >= 90]
                candidate_head_contigs_id = list(head_df[1].values)
                candidate_head_contigs = [contig_dic[x] for x in candidate_head_contigs_id]
                # valueable_contigs.extend(candidate_head_contigs)
                best_head_contig = None
                best_head_contig_score = 0
                for id in candidate_head_contigs_id:
                    head_contig = contig_dic[id]
                    score = findSupportReadScore(head_contig, sequences_scores)
                    if score > best_head_contig_score:
                        best_head_contig = head_contig
                        best_head_contig_score = score
                if best_head_contig not in best_contigs:
                    best_contigs.append(best_head_contig)
                if best_head_contig not in template.best_fragments:
                    template.best_fragments.append(best_head_contig)
                if fragment not in template.best_fragments:
                    template.best_fragments.append(fragment)
                # else:
                #     template.best_fragments.append(fragment)
            except Exception as e:
                print(e)
                quit()

            try:
                tail_df = pd.read_csv(tail_out, delimiter='\t', header=None)
                tail_df = tail_df[tail_df[2] >= 90]
                candidate_tail_contigs_id = list(tail_df[1].values)
                candidate_tail_contigs = [contig_dic[x] for x in candidate_tail_contigs_id]
                # valueable_contigs.extend(candidate_tail_contigs)
                best_tail_contig = None
                best_tail_contig_score = 0
                for id in candidate_tail_contigs_id:
                    tail_contig = contig_dic[id]
                    score = findSupportReadScore(tail_contig, sequences_scores)
                    if score > best_tail_contig_score:
                        best_tail_contig = tail_contig
                        best_tail_contig_score = score
                if best_tail_contig not in best_contigs:
                    best_contigs.append(best_tail_contig)
                if fragment not in template.best_fragments:
                    template.best_fragments.append(fragment)
                if best_tail_contig not in template.best_fragments:
                    template.best_fragments.append(best_tail_contig)
                # hook = best_tail_contig[len(best_tail_contig) - 3:]
                # print(best_tail_contig,hook)
                # hook_out = f'{froot}/hook_refactor.m8'
                # with open(f'{froot}/hook.fasta', 'w') as f:
                #     f.write(f'>tail_hook\n{hook}')
                # os.system(f'rapsearch -q {froot}/hook.fasta -d {froot}/contigs -o {froot}/{froot}_hook')
                # os.system(
                #     f'python processRapsearchM8.py -input {froot}/{froot}_hook.m8 -output {hook_out}')
                # hook_df = pd.read_csv(hook_out, delimiter='\t', header=None)
                #
                # print(hook_df)
            except Exception as e:
                print(e)
                quit()

        # best_result_coverage_list = []
        # is_continue = False
        # for best_result_position in range(len(best_result)):
        #     current = best_result[best_result_position]
        #     if not is_continue and current != ' ':
        #         start = best_result_position
        #         is_continue = True
        #     elif is_continue and current == ' ':
        #         end = best_result_position-1
        #         is_continue = False
        #         best_result_coverage_list.append([start,end])
        #     elif is_continue and len(best_result) == best_result_position + 1:
        #         best_result_coverage_list.append([start,best_result_position])
        # print(best_result_coverage_list)
        template.best_fragments = sorted(template.best_fragments, key=len)
        template.best_fragments = [j for i, j in enumerate(template.best_fragments) if all(j == k or (j not in k) for k in template.best_fragments[i + 1:])]

        step = 250
        html += '*' * 100 + 'Merged Result' + '*' * 100 + '<br>'
        html += 'Template ID: {}<br>'.format(template.id)
        for i in range(0, len(template.sequence), step):
            try:
                sub_template = template.sequence[i:i + step]
            except:
                sub_template = template.sequence[i:]
            # print(sub_template)
            html += '<pre>' + sub_template + '</pre>'
            for sequence in merged_result:
                try:
                    sub_sequence = sequence[i:i + step]
                except:
                    sub_sequence = sequence[i:]
                sub_sequence = ''.join(sub_sequence)
                if not any(c.isalpha() for c in sub_sequence):
                    continue
                html += '<pre>' + sub_sequence + '</pre>'
            html += '<br>'
        for best_contig in best_contigs:
            print(12341234, best_contig)
            # html += '<pre>' + best_contig + '</pre>'
        for best_fragment in template.best_fragments:
            html += '<pre>' + best_fragment + '</pre>'
        html += '<br>'
        html += 'Minimum Contigs Array (Blue part): ' + '<br>'
        for index in range(len(minimum_contigs_array)):
            contig = minimum_contigs_array[index]
            id = uuid.uuid4()
            colored_contig = list(contig.sequence)
            for i in range(contig.contig_interval[0] - 1, contig.contig_interval[1]):
                colored_contig[i] = '<font color="blue">{}</font>'.format(colored_contig[i])
            html += '<pre>' + 'Template interval: ' + str(contig.template_interval) + ' | ' + 'Contig score: ' + str(
                round(findSupportReadScore(contig.sequence, sequences_scores), 2)) + '</pre>'
            html += '<pre>' + ''.join(colored_contig) + '</pre>'
            y = list(contig.rates.values())
            y = list(NormalizeData(y))
            x = [letter for letter in contig.sequence]
            line_chart = f'''
            <canvas id="{id}" style="width:100%;max-width:600px"></canvas>
            <script>
            var xValues = {str(x)};
            var yValues = {str(y)};
            
            new Chart("{id}", 
            '''

            line_chart += '''{
              type: "line",
              data: {
                labels: xValues,
                datasets: [{
                  fill: false,
                  lineTension: 0,
                  backgroundColor: "rgba(0,0,255,1.0)",
                  borderColor: "rgba(0,0,255,0.1)",
                  data: yValues
                }]
              },
              options: {
                legend: {display: false},
              }
            });
            </script>                       
            '''
            html += line_chart + '<br>'

        html += 'unused Reads (Green part): ' + '<br>'
        for read in unused_reads_intervals.keys():
            count = 0
            for i in range(unused_reads_intervals[read][0][0] - 1, unused_reads_intervals[read][0][1]):
                if 'green' in merged_result[0][i]:
                    count += 1
            if count == 0:
                continue
            read_sequence = list(read)
            for i in range(unused_reads_intervals[read][1][0] - 1, unused_reads_intervals[read][1][1]):
                read_sequence[i] = '<font color="green">{}</font>'.format(read_sequence[i])
            html += '<pre>' + 'Template interval: ' + str(
                unused_reads_intervals[read][0]) + ' | ' + 'Read Count: ' + str(reads_count[read]) + ' | ' + ''.join(
                read_sequence) + '</pre>'
            # html += '<pre>' + ''.join(read_sequence) + '</pre>'

    # valueable_contigs = list(Counter(valueable_contigs).keys())
    # graph = naive_db.construct_naive_debruijn_graph(valueable_contigs,4,False)
    # outputs = naive_db.output_contigs(graph,[],[])
    # outputs = sorted(outputs,key=lambda x:findSupportReadScore(x,sequences_scores),reverse=True)
    # for output in outputs:
    #     print(output)
    light = ['','']
    heavy = ['','']
    for Template in Templates:
        if 'Light' in Template.id:
            if 'NonConstant' in Template.id:
                light[0] = Template
            else:
                light[1] = Template
        else:
            if 'NonConstant' in Template.id:
                heavy[0] = Template
            else:
                heavy[1] = Template
    # for chain in light:
    #     inputs = chain.best_fragments
    #     print(inputs)
    #     start_inputs = inputs[0]
    #     start_inputs.extend(inputs[1])
    #     print(start_inputs)
    #
    #     inputs.remove(inputs[0])
    #     inputs.remove(inputs[1])
    #     graph = naive_db.construct_naive_debruijn_graph(start_inputs, 5, False)
    #     outputs = naive_db.output_contigs(graph)
    #     outputs = sorted(outputs, key=lambda x: findSupportReadScore(x, sequences_scores), reverse=True)
    #     for output in outputs:
    #         print(output)
    #     quit()

    html += '''
    </body>
    </html>
    '''
    htmlFile.write(html)
    htmlFile.close()
    outFile.write(message)
    outFile.close()
