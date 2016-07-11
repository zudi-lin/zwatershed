# convolutional network metric scripts
- Code for fast watersheds. Code is based around code from https://bitbucket.org/poozh/watershed described in http://arxiv.org/abs/1505.00249.  For use in https://github.com/naibaf7/PyGreentea. 

# building 

### conda
- `conda install -c conda-forge zwatershed`

### pip [<img src="https://img.shields.io/pypi/v/zwatershed.svg?maxAge=2592000">](https://pypi.python.org/pypi/zwatershed/)
- `pip install zwatershed`
- $3^5$

### from source
- clone the repository
- run ./make.sh

### requirements
- numpy, h5py, cython
- if using parallel watershed, also requires multiprocessing or pyspark
- in order to build the cython, requires a c++ compiler and boost

# function api
- `(segs, rand) = zwatershed_and_metrics(segTrue, aff_graph, eval_thresh_list, seg_save_thresh_list)`
	- *returns segmentations and metrics*
	- `segs`: list of segmentations
		- `len(segs) == len(seg_save_thresh_list)`
	- `rand`: dict
		- `rand['V_Rand']`:  V_Rand score (scalar)
		- `rand['V_Rand_split']`: list of score values
			- `len(rand['V_Rand_split']) == len(eval_thresh_list)`
		- `rand['V_Rand_merge']`: list of score values, 
			- `len(rand['V_Rand_merge']) == len(eval_thresh_list)`
- `segs = zwatershed(aff_graph, seg_save_thresh_list)` 
		- *returns segmentations*
	- `segs`: list of segmentations
		- `len(segs) == len(seg_save_thresh_list)`

##### These methods have versions which save the segmentations to hdf5 files instead of returning them
- `rand = zwatershed_and_metrics_h5(segTrue, aff_graph, eval_thresh_list, seg_save_thresh_list, seg_save_path)`
- `zwatershed_h5(aff_graph, eval_thresh_list, seg_save_path)`

##### All 4 methods have versions which take an edgelist representation of the affinity graph
- `(segs, rand) = zwatershed_and_metrics_arb(segTrue, node1, node2, edgeWeight, eval_thresh_list, seg_save_thresh_list)`
- `segs = zwatershed_arb(seg_shape, node1, node2, edgeWeight, seg_save_thresh_list)`
- `rand = zwatershed_and_metrics_h5_arb(segTrue, node1, node2, edgeWeight, eval_thresh_list, seg_save_thresh_list, seg_save_path)`
- `zwatershed_h5_arb(seg_shape, node1, node2, edgeWeight, eval_thresh_list, seg_save_path)`

# parallel watershed - 4 steps
- *a full example is given in par_ex.ipynb*

1. Partition the subvolumes
	- `partition_data = partition_subvols(pred_file,out_folder,max_len)`
		- evenly divides the data in *pred_file* with the constraint that no dimension of any subvolume is longer than max_len
2. Zwatershed the subvolumes
	1. `eval_with_spark(partition_data[0])`
		- *with spark*
	2. `eval_with_par_map(partition_data[0],NUM_WORKERS)`
		- *with python multiprocessing map*
	- after evaluating, subvolumes will be saved into the out\_folder directory named based on their smallest indices in each dimension (ex. path/to/out\_folder/0\_0\_0\_vol)
3. Stitch the subvolumes together
	- `stitch_and_save(partition_data,outname)`
		- stitch together the subvolumes in partition_data
		- save to the hdf5 file outname
			- outname['starts'] = list of min_indices of each subvolume
			- outname['ends'] = list of max_indices of each subvolume
			- outname['seg'] = full stitched segmentation
			- outname['seg_sizes'] = array of size of each segmentation
			- outname['rg_i'] = region graph for ith subvolume
4. Threshold individual subvolumes by merging
	- `seg_merged = merge_by_thresh(seg,seg_sizes,rg,thresh)`
		- load in these areguments from outname
