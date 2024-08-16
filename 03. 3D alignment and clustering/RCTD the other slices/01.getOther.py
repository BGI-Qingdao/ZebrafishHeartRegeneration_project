import pandas as pd
import anndata as ad

data = ad.read_h5ad('zebrafish.3D.h5ad')
meta = pd.read_csv('scvi_reduction.nice.meta.txt',sep='\t',header=0,index_col=0)
data.obs['nice1.2'] = meta['RNA_snn_res.1.2']

data0 = data[data.obs['nice1.2'].isnull()==False]
nice_sids = data0.obs['sid'].unique()
data1 = data[~data.obs['sid'].isin(nice_sids)]

for sid in data1.obs['sid'].unique():
    sdata = data1[data1.obs['sid']==sid]
    sdata.write(f'other.{sid}.h5ad',compression='gzip')
