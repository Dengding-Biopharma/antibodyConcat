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
    chain = 'Heavy'

    templates = read_fasta('homo_templates.fasta')

    file = open(f'NonConstant_{species}_{chain}.fasta','w')

    for template_key in templates.keys():
        try:
            ann = annotation[template_key]
            template_sequence = templates[template_key]
            print(template_sequence)
            for ann_key in ann.keys():
                if chain in template_key and chain in template_key:
                    if 'CONSTANT' in ann_key:
                        file.write(f'>{template_key}\n{template_sequence[:ann[ann_key][0]]}\n')

        except:
            pass

