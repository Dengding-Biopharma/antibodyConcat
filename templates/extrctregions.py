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
            value = [int(value[0])-1,int(value[1])-1]
            ann[id][key] = value
    return ann

def read_fasta(path,species=None):
    f = open(path, 'r')
    lines = f.readlines()
    f.close()
    dic = {}
    if species:
        for i in range(len(lines)):
            line = lines[i]
            if line[0] == '>' and species in line:
                id = line.rstrip()[1:]
                contig = lines[i + 1]
                dic[id] = contig.rstrip()
        return dic
    else:
        for i in range(len(lines)):
            line = lines[i]
            if line[0] == '>':
                id = line.rstrip()[1:]
                contig = lines[i + 1]
                dic[id] = contig.rstrip()
        return dic

if __name__ == '__main__':
    annotation = read_ann('mAB_database.ann')
    species = 'Homo'
    chain = 'Light'

    templates = read_fasta('homo_templates.fasta')

    fr1_file = open(f'fr1_{species}_{chain}.fasta','w')
    fr2_file = open(f'fr2_{species}_{chain}.fasta','w')
    fr3_file = open(f'fr3_{species}_{chain}.fasta','w')
    fr4_file = open(f'fr4_{species}_{chain}.fasta','w')
    constant_file = open(f'constant_{species}_{chain}.fasta','w')
    for template_key in templates.keys():
        try:
            ann = annotation[template_key]
            template_sequence = templates[template_key]
            print(template_sequence)
            for ann_key in ann.keys():
                if chain in template_key and chain in template_key:
                    if 'FR1' in ann_key:
                        print(ann_key,template_sequence[ann[ann_key][0]:ann[ann_key][1]])
                        fr1_file.write(f'>{template_key}\n{template_sequence[ann[ann_key][0]:ann[ann_key][1]]}\n')
                    if 'FR2' in ann_key:
                        print(ann_key,template_sequence[ann[ann_key][0]:ann[ann_key][1]])
                        fr2_file.write(f'>{template_key}\n{template_sequence[ann[ann_key][0]:ann[ann_key][1]]}\n')
                    if 'FR3' in ann_key:
                        print(ann_key,template_sequence[ann[ann_key][0]:ann[ann_key][1]])
                        fr3_file.write(f'>{template_key}\n{template_sequence[ann[ann_key][0]:ann[ann_key][1]]}\n')
                    if 'FR4' in ann_key:
                        print(ann_key,template_sequence[ann[ann_key][0]:ann[ann_key][1]])
                        fr4_file.write(f'>{template_key}\n{template_sequence[ann[ann_key][0]:ann[ann_key][1]]}\n')
                    if 'CONSTANT' in ann_key:
                        print(ann_key, template_sequence[ann[ann_key][0]:ann[ann_key][1]])
                        constant_file.write(f'>{template_key}\n{template_sequence[ann[ann_key][0]:ann[ann_key][1]]}\n')
        except:
            pass


