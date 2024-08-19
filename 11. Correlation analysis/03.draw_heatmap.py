import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import numpy as np
plt.rcParams['font.size'] = 15
infile = sys.argv[1]
prefix = sys.argv[2]

data = pd.read_csv(infile,sep='\t',header=0,index_col=0)
data = data.replace(np.inf,0)
data = data.replace(np.nan,0)
data = data.T

from matplotlib.colors import LinearSegmentedColormap
colorfull = LinearSegmentedColormap.from_list('colorful',
                ['#adb5c5','#3e61ab','#6bc9e8','#f5e80b','#ff8000','#e3080b'],N=256)
sns.clustermap(data,figsize=(10,13.5), vmin=0, cmap='binary')
plt.savefig(f'{prefix}.pdf')
