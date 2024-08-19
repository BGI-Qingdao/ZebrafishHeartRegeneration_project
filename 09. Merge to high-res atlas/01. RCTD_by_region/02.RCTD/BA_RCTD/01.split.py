import anndata as ad
import sys

adata = ad.read('BA.h5ad')
for x in adata.obs['sid'].unique():
    tmp = adata[ adata.obs['sid']==x ]
    tmp.write(f'00.data/{x}.h5ad',compression='gzip')
