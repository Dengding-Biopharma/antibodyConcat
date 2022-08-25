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


def toHTML(string):
    result = []
    for char in string:
        if char == ' ':
            result.append('-')
        else:
            result.append(char)
    result = ''.join(result)

    return result
def get_args():
    parser = argparse.ArgumentParser()
    # start
    parser.add_argument('-froot', type=str,required=True)
    parser.add_argument('-template', type=str,required=True)
    parser.add_argument('-source', type=str,required=True)
    args = parser.parse_args()
    return args

def NormalizeData(data):
    return (data - np.min(data)) / (np.max(data) - np.min(data))


if __name__ == '__main__':
    args = get_args()
    template = args.template
    template_name = f'templates/{template}_templates.fasta'
    froot = args.froot
    contig_filepath = f'{froot}/{froot}_sorted.fasta'

    settingFile = open(f'{froot}/setting.json', 'r')
    setting = json.load(settingFile)
    sourceFilePath = args.source
    DF = pd.DataFrame()
    unused_reads = pd.DataFrame()
    sequences_scores = dict()
    for root, dir, files in os.walk(sourceFilePath):
        root = root + '/'
        for file in files:
            filename = root + file
            data = pd.read_csv(filename, delimiter='\t')

            unused_reads_temp = data[data['Score'] >= 0.1]
            unused_reads_temp = unused_reads_temp[-50 < unused_reads_temp['PPM Difference']]
            unused_reads_temp = unused_reads_temp[unused_reads_temp['PPM Difference'] < 50]
            unused_reads_temp.reset_index(inplace=True, drop=True)
            unused_reads = unused_reads.append(unused_reads_temp)

            temp = data[data['Score'] >= setting['score_cut']]
            temp = temp[-50 < temp['PPM Difference']]
            temp = temp[temp['PPM Difference'] < 50]
            temp.reset_index(inplace=True, drop=True)
            for i in range(len(temp)):
                if temp['DENOVO'][i] not in sequences_scores.keys():
                    sequences_scores[temp['DENOVO'][i]] = temp['Score'][i]
                else:
                    sequences_scores[temp['DENOVO'][i]] = temp['Score'][i] + sequences_scores[temp['DENOVO'][i]]
            DF = DF.append(temp)

    DF.reset_index(inplace=True, drop=True)
    unused_reads.reset_index(inplace=True,drop=True)

    os.system(f'prerapsearch -d {template_name} -n templates/temp-db')
    os.system(f'rapsearch -q {froot}/{froot}_sorted.fasta -d templates/temp-db -o {froot}/rapsearch_outputs -z 6')
    os.system(f'python processRapsearchM8.py -input {froot}/rapsearch_outputs.m8 -output {froot}/rapsearch_outputs_refactor.m8')
    df = pd.read_csv(f'{froot}/rapsearch_outputs_refactor.m8',delimiter='\t',header=None)
    # df = pd.read_csv(f'{froot}/{froot}_blasthomoTemplate.m8',delimiter='\t',header=None)
    df = df[df[2] >= 80]
    df = df.sort_values(by=0)
    df = df.reset_index(drop=True)

    template_dic = read_fasta(template_name)
    templates = list(template_dic.keys())
    contig_dic = read_fasta(contig_filepath)
    contigs = list(contig_dic.keys())

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

    for template_id in template_contig_group.keys():
        template = Template(template_id, template_dic[template_id])
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
            if template.contigArrays[i][0].template_interval[0]<start_contig.template_interval[0] \
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
                    if current_contig.template_interval[0] < contig.template_interval[0] < current_contig.template_interval[1] < \
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
            unused_reads_intervals[unusedReads_dic[unusedRead]] = [[template_left,template_right],[read_left,read_right]]
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

        coverage = matched_length/len(template.sequence)
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
        html += 'Minimum Contigs Array (Blue part): ' + '<br>'
        for index in range(len(minimum_contigs_array)):
            contig = minimum_contigs_array[index]
            id = uuid.uuid4()
            colored_contig = list(contig.sequence)
            for i in range(contig.contig_interval[0] - 1,contig.contig_interval[1]):
                colored_contig[i] = '<font color="blue">{}</font>'.format(colored_contig[i])
            html += '<pre>' + 'Template interval: '+ str(contig.template_interval) + ' | ' + 'Contig score: ' + str(round(findSupportReadScore(contig.sequence,sequences_scores),2)) + '</pre>'
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
            for i in range(unused_reads_intervals[read][0][0]-1,unused_reads_intervals[read][0][1]):
                if 'green' in merged_result[0][i]:
                    count += 1
            if count == 0:
                continue
            read_sequence = list(read)
            for i in range(unused_reads_intervals[read][1][0]-1,unused_reads_intervals[read][1][1]):
                read_sequence[i] = '<font color="green">{}</font>'.format(read_sequence[i])
            html += '<pre>' + 'Template interval: ' + str(unused_reads_intervals[read][0]) + ' | ' + 'Read Count: ' + str(reads_count[read]) + ' | ' + ''.join(read_sequence) + '</pre>'
            # html += '<pre>' + ''.join(read_sequence) + '</pre>'

    html += '''
    </body>
    </html>
    '''
    htmlFile.write(html)
    htmlFile.close()
    outFile.write(message)
    outFile.close()
