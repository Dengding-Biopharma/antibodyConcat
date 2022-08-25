import argparse

import pandas as pd
from tqdm import trange


def get_args():
    parser = argparse.ArgumentParser()
    # start
    parser.add_argument('-input', type=str)
    parser.add_argument('-output',type=str)
    args = parser.parse_args()
    return args

args = get_args()
path = args.input

file = open(path,mode='r+')
lines = file.readlines()
file.close()
df = pd.DataFrame(columns=[0,1,2,3,4,5,6,7,8,9,10,11])
outlines = []
aline = []
for i in trange(len(lines)):
    line = lines[i]
    if '#' not in line:
        items = line.split('\t')
        for item in items:
            if item:
                aline.append(item.rstrip())
    if len(aline) == 12:
        outlines.append(aline)
        aline = []

df = pd.DataFrame(outlines,columns=[0,1,2,3,4,5,6,7,8,9,10,11])
out = args.output
df.to_csv(out,index=False,sep='\t',header=False)