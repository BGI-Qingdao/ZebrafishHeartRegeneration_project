import sys
import re
import os
import pandas as pd
import numpy as np
import scvelo as scv

scv.logging.print_version()
scv.settings.verbosity = 3
scv.settings.presenter_view = True

loom = scv.read(sys.argv[1],cache=False)
loom.var_names_make_unique()
loom.obs = loom.obs.rename(index = lambda x: re.sub(":","_",x))
loom.obs.head()

outpre = sys.argv[3]

meta = pd.read_table(sys.argv[2])
select=meta

adata = loom[np.isin(loom.obs.index, select["CellID"])]
select_merge = pd.DataFrame(adata.obs.index).merge(select,on="CellID")
adata.obsm['X_umap'] = select_merge[["umap_1","umap_2"]].values
select_merge = select_merge.set_index('CellID')
adata.obs['celltype'] = select_merge["celltype"]
adata.uns['clusters'] = select_merge[["celltype"]].values

adata.var_names_make_unique()

scv.pp.filter_and_normalize(adata,min_shared_counts=30, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)
adata.write(outpre + '.dynamical.h5ad', compression = 'gzip')

ident_colours = {
"#E31A1CFF","#33A02CFF","#006400","#C0FF3E","#B2DF8AFF","#d2fad2","#91D0BE","#68A180","#C5DEBA","#FDBF6FFF",
"#FF7F00FF","#A6CEE3FF","#1F78B4FF","#57C3F3","#FB9A99FF","#E95C59","#B53E2B","#014d01","#03ad03","#E5D2DD",
"#6A3D9AFF","#E59CC4","#CAB2D6FF","#DA70D6","#E4C755","#FFD700","#fff8d1","#FFFF99FF","#D2B48C","#F5DEB3",
"#F5DEB3","#8B4726","#FFFACD","#EEE9BF","#E8E8E8","#F0FFF0","#000000","#363636","#CDC9C9","#8B8989","#4d4c4c"
}

scv.pl.velocity_embedding_stream(adata, save= outpre + ".scv1.pdf", basis='X_umap',color = "celltype", palette = ident_colours, alpha = 1, legend_loc="right margin", size=10)
scv.pl.velocity_embedding_stream(adata, save= outpre + ".scv2.pdf", basis='X_umap',color = "celltype", palette = ident_colours, alpha = 1, legend_loc="on data", size=10, legend_fontsize=5)

