
def read_fasta(fname):
    f = open(fname, 'r')
    lines = f.readlines()
    f.close()
    fasta={}

    for i in range(len(lines)):
        line = lines[i]
        if '9844' in line:
            fasta[line] = lines[i+1]
    return fasta

fasta = read_fasta('nanobody_template.fasta')


outfile = 'nanobody_templates.fasta'
f = open(outfile,'w')
for key in fasta.keys():
    f.writelines(key)
    f.writelines(fasta[key])
f.close()