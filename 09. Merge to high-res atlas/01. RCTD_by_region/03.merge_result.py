import pandas as pd
import sys
import os
alldata = []
for i in range(200):
    tmp = []
    for f in ['A','BA','V']:
        if not os.path.isfile(f'{f}/s{i}.result.csv'):
            continue
        tmp.append(pd.read_csv(f'{f}/s{i}.result.csv',sep=',',header=0))
    if len(tmp)>0:
        tmp = pd.concat(tmp,ignore_index=True)
        tmp.to_csv('s{i}.result.csv')