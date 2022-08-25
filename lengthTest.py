from collections import Counter

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt


file_6h = pd.read_csv('test/TS-HCl-P-6h.denovo.csv', sep='\t')
file_30min = pd.read_csv('test/TS-HCl-P-30min.denovo.csv', sep='\t')

file_6h = file_6h[file_6h['Score'] >= 0.8]
file_30min = file_30min[file_30min['Score'] >= 0.8]
denovo_6h = list(Counter(file_6h['DENOVO'].values))
denovo_30min = list(Counter(file_30min['DENOVO'].values))

venn2([set(denovo_6h),set(denovo_30min)],set_labels=('P-6h','P-30min'))
plt.title('TS-HCL-P')
plt.show()
common_set = set(denovo_6h) & set(denovo_30min)
denovo_6h_unique = list(set(denovo_6h) - common_set)
denovo_30min_unique = list(set(denovo_30min) - common_set)
denovo_6h_length = [len(x) for x in denovo_6h_unique]
denovo_30min_length = [len(x) for x in denovo_30min_unique]
denovo_6h_length_dic = dict()
for length in denovo_6h_length:
    if length not in denovo_6h_length_dic.keys():
        denovo_6h_length_dic[length] = 1
    else:
        denovo_6h_length_dic[length] += 1
denovo_6h_length_dic = dict(sorted(denovo_6h_length_dic.items(),key=lambda x:x[0]))
denovo_30min_length_dic = dict()
for length in denovo_30min_length:
    if length not in denovo_30min_length_dic.keys():
        denovo_30min_length_dic[length] = 1
    else:
        denovo_30min_length_dic[length] += 1
denovo_30min_length_dic = dict(sorted(denovo_30min_length_dic.items(),key=lambda x:x[0]))
print(denovo_30min_length_dic)
print(denovo_6h_length_dic)
maxlength = max(max(list(denovo_30min_length_dic.keys()),list(denovo_6h_length_dic.keys())))
print(maxlength)

bins = np.linspace(1, maxlength,maxlength)
print(bins)

x = []
for length in range(1,maxlength+1):
    if length in denovo_30min_length_dic.keys():
        x.append(denovo_30min_length_dic[length])
    else:
        x.append(0)
x = np.array(x)
y = []
for length in range(1,maxlength+1):
    if length in denovo_6h_length_dic.keys():
        y.append(denovo_6h_length_dic[length])
    else:
        y.append(0)
y = np.array(y)
X_axis = np.arange(len(bins))
plt.bar(X_axis - 0.2,x,0.4,label='30min')
plt.bar(X_axis + 0.2,y,0.4,label='6h')

plt.xticks(X_axis,bins)
plt.title('TS-HCL-P')
plt.xlabel("Length")
plt.ylabel("Number of counts")
plt.legend()
plt.show()
