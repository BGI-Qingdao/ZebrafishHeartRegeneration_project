import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import sys

in1 = sys.argv[1]
in2 = sys.argv[2]
prefix = sys.argv[3]

data1 = pd.read_table(in1, sep='\t',header=0, index_col=0)
#data1 = data1.replace(np.inf, 0)
#data1 = data1.replace(np.nan, 0)
data2 = pd.read_table(in2, sep='\t',header=0, index_col=0)
#data2 = data2.replace(np.inf, 0)
#data2 = data2.replace(np.nan, 0)

common_gene = np.intersect1d(data1.index.to_numpy(),data2.index.to_numpy())
print(f'#gene from {in1}: {len(data1.index)}')
print(f'#gene from {in2}: {len(data2.index)}')
print(f'#common gene : {len(common_gene)}',flush=True)
data1 = data1.loc[common_gene].copy()
data2 = data2.loc[common_gene].copy()

category1 = data1.columns.tolist()
category2 = data2.columns.tolist()

sims =  np.zeros((len(category2),len(category1)),dtype=float)

for i,c2 in enumerate(category2):
    for j,c1 in enumerate(category1):
        #sims[i,j],_ =  pearsonr(data1[c1].to_numpy(),data2[c2].to_numpy())
        tmp = np.corrcoef(data1[c1].to_numpy(),data2[c2].to_numpy())
        sims[i,j] = tmp[0,1]
ret = pd.DataFrame(data=sims,columns=category1,index=category2)
ret.to_csv(f'{prefix}.heatmap.csv',header=True,index=True,sep='\t')
