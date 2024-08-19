import pandas as pd
import anndata as ad
import numpy as np
import sys,getopt
import scipy

def grouped_obs_mean(adata, group_key, layer=None, gene_symbols=None):
    if layer is not None:
        getX = lambda x: x.layers[layer]
    else:
        getX = lambda x: x.X

    if gene_symbols is not None:
        new_idx = adata.var[idx]
    else:
        new_idx = adata.var_names
    grouped = adata.obs.groupby(group_key)
    out = pd.DataFrame(
        np.zeros((adata.shape[1], len(grouped)), dtype=np.float64),
        columns=list(grouped.groups.keys()),
        index=adata.var_names
    )

    for group, idx in grouped.indices.items():
        X = getX(adata[idx])
        if scipy.sparse.issparse(X):
            X = X.todense()
        out[group] = np.array(np.ravel(X.mean(axis=0, dtype=np.float64))).tolist()
    return out

def usage():
    print("""
python3 average.py -i <xxx.input.h5ad>
                   -k <anno_key>
                   -o <prefix>
                   -n [optional, norm total or not, default not]
                   -p [optional, log1p  or not, default not]
                   -s [optional, scale or not, default not]
""")

def main(argv):
    input_h5ad = ''
    prefix = ''
    anno = ''
    norm_it = False
    scale_it = False
    log_it = False
    try:
        opts, args = getopt.getopt(argv,"hi:o:k:nsp",["help",])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h' ,'--help'):
            usage()
            sys.exit(0)
        elif opt in ("-o", ):
            prefix = arg
        elif opt in ("-i", ):
            input_h5ad = arg
        elif opt in ("-k"):
            anno = arg
        elif opt in ("-n"):
            norm_it = True
        elif opt in ("-p"):
            log_it = True
        elif opt in ("-s"):
            scale_it = True

    if input_h5ad == '' or anno == '' or prefix == '':
        usage()
        sys.exit(0)

    adata = ad.read_h5ad(input_h5ad)
    if norm_it:
        import scanpy as sc
        print('now call sc.pp.normalize_total',flush=True)
        sc.pp.normalize_total(adata)
    if log_it:
        import scanpy as sc
        print('now call sc.pp.log1p',flush=True)
        sc.pp.log1p(adata)
    if scale_it:
        import scanpy as sc
        print('now call sc.pp.scale',flush=True)
        sc.pp.scale(adata, max_value=6)
    ret = grouped_obs_mean(adata,anno)
    ret.to_csv(f'{prefix}.average.csv',sep='\t',header=True,index=True)

if __name__ == '__main__':
    main(sys.argv[1:])
