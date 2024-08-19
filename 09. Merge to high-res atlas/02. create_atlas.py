import anndata as ad
import numpy as np
import pandas as  pd
import sys
import os
########################################
# load RCTD result
########################################
alldata = []
for i in range(200):
    fname = f'00.RCTD/s{i}.result.csv'
    if os.path.exists(fname):
        tmp = pd.read_csv(fname,sep=',',header=0,index_col=0,quotechar='"')
        tmp['cell'] = tmp.index.to_list()
        tmp['cell'] = tmp.apply(lambda row: f's{i}_{row["cell"]}',axis=1)
        tmp = tmp[['cell',"spot_class","first_type","second_type"]].copy()
        tmp = tmp[tmp['spot_class'].isin(['singlet','doublet_certain'])].copy()
        alldata.append(tmp)

alldata = pd.concat(alldata,ignore_index=True)
alldata.to_csv('all.csv',sep='\t',header=True,index=False)
alldata = alldata.set_index('cell')
adata = ad.read_h5ad('ST.3D.count.region.h5ad')
alldata['region'] = adata.obs['region']
alldata['3D_x'] = adata.obs['3D_x']
alldata['3D_y'] = adata.obs['3D_y']
alldata['3D_z'] = adata.obs['3D_z']
alldata.to_csv('all.csv',sep='\t',header=True,index=True)
########################################
# create pseudo spot 
########################################
ctdata = alldata[["first_type",'region','3D_x','3D_y','3D_z']].copy()
ctdata.columns = ['celltype','region','3D_x','3D_y','3D_z']
dump = alldata[alldata['spot_class']=='doublet_certain'].copy()
dump = dump[["second_type",'region','3D_x','3D_y','3D_z']].copy()
dump.columns = ['celltype','region','3D_x','3D_y','3D_z']
dump['3D_x'] = dump['3D_x']+10
dump['3D_z'] = dump['3D_z']+10

ctdata['cell'] = ctdata.index
dump['cell'] = dump.index

newdata = pd.concat([ctdata,dump],ignore_index=True)
newdata.to_csv('atlas.csv',sep='\t',header=True,index=False)

########################################
# merge scRNA & ST
########################################
dump = pd.read_csv('atlas.csv',sep='\t',header=0)
dump=dump[dump['celltype']!="Cardiomyocytes (A)"]
dump=dump[dump['celltype']!="Cardiomyocytes (V)"]
dump=dump[dump['celltype']!="Red blood cells"]
dump = dump[['celltype',	'region','3D_x','3D_y','3D_z']].copy()

stdata = pd.read_csv('ST.atlas.csv',sep='\t',header=0)
stdata['cell_type_all'].unique().tolist()
stdata=stdata[stdata['cell_type_all'].isin(['Epicardium','Red blood cells','Cardiomyocytes (Atrium)','Endothelial cells','Cardiomyocytes (Ventricle)', 'Valves','Nerve cells'])]
stdata = stdata[['cell_type_all','region','3D_x','3D_y','3D_z']].copy()
stdata.columns = ['celltype',	'region','3D_x','3D_y','3D_z']

alldata = pd.concat([stdata,dump],ignore_index=True)

########################################
# save atlas h5ad
########################################
ngene = 2
fakeX = np.zeros((len(alldata),ngene),dtype=float)
adata = ad.AnnData(X=fakeX,obs=alldata)


adata.obsm['spatial'] = adata.obs[['3D_x','3D_y','3D_z']].to_numpy()
adata.write('atlas.h5ad',compression='gzip')