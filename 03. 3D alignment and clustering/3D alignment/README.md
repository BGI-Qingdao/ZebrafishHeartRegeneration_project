## The pipeline of 3D alignment

* step1: alignment the ssDNA images of slices via TrakEM2
* output1: trakEM2.json 

* step2: prepare config file for step3

```bash
#!/bin/bash

function get_xmin(){
    local infile=$1
    grep final $infile|awk -F '=|,| ' '{print $6}'
}

function get_ymin(){
    local infile=$1
    grep final $infile|awk -F '=|,| ' '{print $10}'
}

function get_affine(){
    local ssid=$1
    temp=""
    if [[ ${#ssid} == "2" ]] ; then
        temp="mask.0"${ssid}".tif"
    else
        temp="mask."${ssid}".tif"
    fi
    grep -A3 $temp trakEM2.json |tail -n1 |awk -F "\"" '{printf("\"%s\"\n",$2);}'
}

H5path='path/zebrafish_heart/04.cluster_bin50/00.h5ad/'
maskpath='path/zebrafish_heart/01.mask_gem/'

echo "flag,h5ad,Z_values,3D_forward,x_shift,y_shift"
for sprefix in `cat bin50.in.csv`
do
    sid=${sprefix:18}
    flag="s"$sid
    h5ad=${H5path}/s${sid}.fu100.h5ad
    Zvalues=$sid
    mask_log=${maskpath}/${sprefix}/mask.log
    xshift=`get_xmin $mask_log`
    yshift=`get_ymin $mask_log`
    D3forward=`get_affine $sid`
    echo $flag","$h5ad","$Zvalues","$D3forward","$xshift","$yshift
done
```
* output1: bin50.3D.csv

* step3: merge slices into one 3D dataset

please download GEM3D_TK from https://github.com/BGI-Qingdao/GEM3D_toolkit

```bash
 python3 GEM3D_toolkit/GEM_toolkit.py apply_alignment -i bin50.the3D.csv -o zebrafish.3D -m True -S True
```
* output3 zebrafish.3D.h5ad