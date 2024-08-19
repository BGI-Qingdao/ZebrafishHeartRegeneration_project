# Run RCTD seperatelly.

* step1, split the stereo-serq 3D altas into A, V, and BA atlas by the mesh model
* step2, run RCTD by corresponding scRNA data 
* step3, merge the results of each slice into one final result.

## step1:

please download VT3D from https://github.com/BGI-Qingdao/VT3D

```bash
python3 vt3d Auxiliary  MeshSub -i zebrafish.3D.h5ad -m atrial1.obj -o A.h5ad --spatial_key spatial
python3 vt3d Auxiliary  MeshSub -i zebrafish.3D.h5ad -m ventricular1.obj -o V.h5ad --spatial_key spatial
python3 vt3d Auxiliary  MeshSub -i zebrafish.3D.h5ad -m arterial_bulb1.obj -o BA.h5ad --spatial_key spatial
```

