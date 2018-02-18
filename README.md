# Cython version of Zwatershed

## Resource
- paper: http://arxiv.org/abs/1505.00249
- original C++ version: https://bitbucket.org/poozh/watershed 
- cython version: https://github.com/jakirkham/zwatershed
- julia version: https://github.com/seung-lab/Watershed.jl
- **ours**: start from the cython version and modify according to the julia
  version (with a bug fix)

## Install
`python setup.py install`

## Package Functions
representation: affinity -> segmentation -> region graph

1. main functions:
    - zwatershed: ```zwatershed()```

2. initial segmentation: ```zw_initial(affs, affs_low, affs_high)```
    - steepest ascent: ```zw_steepest_ascent(aff, low, high)```
    - divide plateaus: ```zw_divide_plateaus(seg)```
    - find basin:  ```zw_find_basins(seg)```

3. merge region:  
    - get region graph: ```zw_get_region_graph(aff, seg) ```
    - merge (with or without mst)
    ```zw_merge_segments_with_function(seg, rg_affs, id1, id2, counts, T_size, T_weight, T_dust, T_merge, T_mst)``` 
    - mean-affinity: 

## Example Scripts
- zwatershed: ```do_zwatershed.py```
- mean-affinity agglomeration: ```mean_aff.py```
