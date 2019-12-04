# Cython version of Zwatershed

## The py3 branch of zwatershed allows you to run the demo notebook (pytorch_connectomics/demo/segmentation.ipynb) with python3
1. A small change has been made to the file zwatershed.pyx by adding the byte encoding tag to dictionary keys
     - fixed error message related to string and bytes encoding in cython, python 2 and python 3:
     - Details about the issue: https://github.com/PySlurm/pyslurm/wiki/Strings-and-bytes-in-Cython
     - example fix:
`cdef np.ndarray[uint64_t, ndim=1] in_seg = np.array(map_init['seg'],dtype='uint64') -->  cdef np.ndarray[uint64_t, ndim=1] in_seg = np.array(map_init[b'seg'],dtype='uint64')`


## Install
`conda install --yes --file requirements.txt` or `pip install -r requirements.txt`

`python setup.py install`

## Package Functions
representation: affinity -> segmentation -> region graph

1. end-to-end functions:
```zwatershed(T_threshes=800, T_aff=[0.01,0.8,0.2], T_aff_relative=True, T_dust=600, T_merge=0)```

2. step-by-step for parameter search

- step 1: initial segmentation: ```zw_initial(affs, affs_low, affs_high)```
    - steepest ascent: ```zw_steepest_ascent(aff, affs_low, affs_high)```
    - divide plateaus: ```zw_divide_plateaus(seg)```
    - find basin:  ```zw_find_basins(seg)```

- step 2: merge region:  
    - get region graph: ```zw_get_region_graph(aff, seg) ```
    - merge (with or without mst)
    ```zw_merge_segments_with_function(seg, rg_affs, id1, id2, counts, T_size, T_weight, T_dust, T_merge, T_mst)``` 
    - mean-affinity: 

## Test Scripts
- zwatershed: ```do_zwatershed.py```
- mean-affinity agglomeration: ```mean_aff.py```

## Resource
- paper: http://arxiv.org/abs/1505.00249
- original C++ version: https://bitbucket.org/poozh/watershed 
- cython version: https://github.com/jakirkham/zwatershed
- julia version: https://github.com/seung-lab/Watershed.jl
- **ours**: start from the cython version and modify according to the julia
  version (with a bug fix)
