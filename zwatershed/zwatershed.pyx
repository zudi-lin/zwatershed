from libcpp.list cimport list
from libcpp.map cimport map
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.pair cimport pair
from libc.stdint cimport uint64_t
import numpy as np
import os
cimport numpy as np
import h5py
from cpython.object cimport PyObject
from cython.operator cimport dereference as deref, preincrement as inc

#-------------- interface methods --------------------------------------------------------------
def zwatershed(np.ndarray[np.float32_t, ndim=4] affs, 
                  T_threshes, T_aff=[0.01,0.8,0.2], T_aff_relative=True, 
                  T_dust=600, T_merge=0.5,
                  save_path=None):
    # aff stats
    affs = np.asfortranarray(np.transpose(affs, (1, 2, 3, 0)))
    dims = affs.shape
    affs_thres = np.percentile(affs, [t*100 for t in T_aff]) if T_aff_relative else T_aff
    print("1. affinity threshold: ", affs_thres)

    print("2. get initial seg")
    map_init = zw_initial_cpp(dims[0], dims[1], dims[2], &affs[0, 0, 0, 0], affs_thres[0], affs_thres[1])

    cdef np.ndarray[uint64_t, ndim=1] in_seg = np.array(map_init['seg'],dtype='uint64')
    cdef np.ndarray[uint64_t, ndim=1] in_counts = np.array(map_init['counts'],dtype='uint64')
    cdef np.ndarray[np.float32_t, ndim=2] in_rg = np.array(map_init['rg'], dtype='float32').reshape(-1, 3)

    # get segs, stats
    T_threshes.sort()
    out_seg = [None]*len(T_threshes)
    out_rg = [None]*len(T_threshes)
    out_counts = [None]*len(T_threshes)
    for i in range(len(T_threshes)):
        print("3. do thres: ", T_threshes[i], T_dust)
        if(len(in_rg) > 0):
            map_init = merge_region(
                dims[0], dims[1], dims[2], &in_rg[0, 0],
                in_rg.shape[0], &in_seg[0], &in_counts[0], 
                len(in_counts), T_threshes[i], affs_thres[2], T_dust, T_merge)
        in_seg = np.array(map_init['seg'], dtype='uint64')
        in_rg = np.array(map_init['rg'], dtype='float32').reshape(-1, 3)
        in_counts = np.array(map_init['counts'], dtype='uint64')

        seg = in_seg.reshape((dims[2], dims[1], dims[0])).transpose(2, 1, 0)
        if save_path is None:
            out_seg[i] = seg.copy()
            out_rg[i] = in_rg.copy()
            out_counts[i] = in_counts.copy()
        else: # in case out of memory
            f = h5py.File(save_path + '_' + str(T_threshes[i]) + '.h5', 'w')
            ds = f.create_dataset('seg', seg.shape, compression="gzip", dtype=np.uint32)
            ds[:] = seg.astype(np.uint32)
            ds2 = f.create_dataset('rg', seg.shape, compression="gzip", dtype=np.uint32)
            ds2[:] = in_rg.astype(np.uint32)
            ds3 = f.create_dataset('counts', seg.shape, compression="gzip", dtype=np.uint32)
            ds3[:] = in_counts.astype(np.uint32)
            f.close()
        print("\t number of regions: "+ str(len(np.unique(seg))))
    if save_path is None:
        return out_seg, out_rg, out_counts

#################################################
# auxilary function for debug purpose

def zw_initial(np.ndarray[np.float32_t, ndim=4] affs, affs_low, affs_high):
    cdef np.ndarray[uint64_t, ndim=1] counts = np.empty(1, dtype='uint64')
    # affs: z*y*x*3
    affs = np.asfortranarray(np.transpose(affs, (1, 2, 3, 0)))
    dims = affs.shape
    map_init = zw_initial_cpp(dims[0], dims[1], dims[2], &affs[0, 0, 0, 0], float(affs_low), float(affs_high))

    graph = np.array(map_init['rg'], dtype='float32')
    rgn_graph = graph.reshape(len(graph) / 3, 3)
    # for output: seg.shape=[z*y*x]
    seg = np.array(map_init['seg'], dtype='uint64').reshape((dims[2], dims[1], dims[0])).transpose(2,1,0)

    return {'rg': rgn_graph, 'seg': seg,
            'counts': np.array(map['counts'], dtype='uint64')}

def zw_merge_region(np.ndarray[uint64_t, ndim=3] seg_init, np.ndarray[uint64_t, ndim=1] counts, 
                   np.ndarray[np.float32_t, ndim=2] rg, T_threshes, T_aff_merge=0.2, 
                    T_dust=600, T_merge=0.5, save_path=None):
    # get initial rg
    # input seg: z*y*x
    # internal seg: x*y*z
    dims = seg_init.shape
    cdef np.ndarray[uint64_t, ndim=1] seg_in = seg_init.transpose(2,1,0).reshape(-1)
    cdef np.ndarray[uint64_t, ndim=1] counts_out = counts
    cdef np.ndarray[np.float32_t, ndim=2] rgn_graph = rg

    # get segs, stats
    T_threshes.sort()
    segs = [None]*len(T_threshes)
    counts_len = len(counts)
    for i in range(len(T_threshes)):
        print "3. do thres: ", T_threshes[i], T_dust
        if(len(rgn_graph) > 0):
            map_init = merge_region(
                dims[0], dims[1], dims[2], &rgn_graph[0, 0],
                rgn_graph.shape[0], &seg_in[0], &counts_out[0], 
                counts_len, T_threshes[i], T_aff_merge, T_dust, T_merge)
            # for next iteration
            seg_in = np.array(map_init['seg'], dtype='uint64')
            counts_out = np.array(map_init['counts'], dtype='uint64')
            counts_len = len(counts_out)
            graph = np.array(map_init['rg'], dtype='float32')
            rgn_graph = graph.reshape(len(graph) / 3, 3)

            seg = seg_in.reshape((dims[2], dims[1], dims[0])).transpose(2, 1, 0)
            if save_path is None:
                segs[i] = seg.copy()
            else:
                f = h5py.File(save_path + '_' + str(T_threshes[i]) + '.h5', 'w')
                ds = f.create_dataset('stack', seg.shape, compression="gzip", dtype=np.uint64)
                ds[:] = seg
                f.close()
            print "\t number of regions: "+ str(len(np.unique(seg)))
        if save_path is None:
            return segs

# step-by-step
def zw_steepest_ascent(np.ndarray[np.float32_t, ndim=4] aff,
                       low, high):
    cdef:
        PyObject *paff = <PyObject *>aff
        PyObject *pseg;
    # input: x*y*z*3
    seg = np.zeros((aff.shape[0], aff.shape[1], aff.shape[2]), np.uint64, order='F')
    pseg = <PyObject *>seg
    steepest_ascent(paff, pseg, low, high);
    return seg

def zw_divide_plateaus(np.ndarray[np.uint64_t, ndim=3] seg):
    cdef:
        PyObject *pseg = <PyObject *>seg
    divide_plateaus(pseg);

def zw_find_basins(np.ndarray[np.uint64_t, ndim=3] seg):
    cdef:
        PyObject *pseg=<PyObject *>seg
        vector[uint64_t] counts
        size_t i
    
    find_basins(pseg, counts)
    print counts.size()
    pycounts = np.zeros(counts.size(), np.uint64)
    for i in range(len(pycounts)):
        pycounts[i] = counts[i]
    return pycounts

def zw_get_region_graph(np.ndarray[np.float32_t, ndim=4] aff,
                        np.ndarray[np.uint64_t, ndim=3] seg):
    '''Return the initial region graph using max edge affinity.
    
    :param aff: the affinity predictions - an array of x, y, z, c where c == 0
                is the affinity prediction for x, c == 1 is the affinity
                prediction for y and c == 2 is the affinity prediction for z
    :param seg: the segmentation after finding basins
    :returns: a region graph as a 3-tuple of numpy 1-d arrays of affinity, 
              ID1 and ID2
    '''
    cdef:
        PyObject *pseg=<PyObject *>seg
        PyObject *paff=<PyObject *>aff
        vector[float] rg_affs
        vector[uint64_t] id1
        vector[uint64_t] id2
        np.ndarray[np.float32_t, ndim=1] np_rg_affs
        np.ndarray[np.uint64_t, ndim=1] np_id1
        np.ndarray[np.uint64_t, ndim=1] np_id2
        size_t i
    max_segid = seg.max()
    get_region_graph(paff, pseg, max_segid, rg_affs, id1, id2)
    np_rg_affs = np.zeros(rg_affs.size(), np.float32)
    np_id1 = np.zeros(id1.size(), np.uint64)
    np_id2 = np.zeros(id2.size(), np.uint64)
    for 0 <= i < rg_affs.size():
        np_rg_affs[i] = rg_affs[i]
        np_id1[i] = id1[i]
        np_id2[i] = id2[i]
    return (np_rg_affs, np_id1, np_id2)

def zw_get_region_graph_average(np.ndarray[np.float32_t, ndim=4] aff,
                                np.ndarray[np.uint64_t, ndim=3] seg):
    '''Return the initial region graph using average edge affinity
    
    :param aff: the affinity predictions - an array of x, y, z, c where c == 0
                is the affinity prediction for x, c == 1 is the affinity
                prediction for y and c == 2 is the affinity prediction for z
    :param seg: the segmentation after finding basins
    :returns: a region graph as a 3-tuple of numpy 1-d arrays of affinity, 
              ID1 and ID2
    '''
    cdef:
        PyObject *pseg=<PyObject *>seg
        PyObject *paff=<PyObject *>aff
        vector[float] rg_affs
        vector[uint64_t] id1
        vector[uint64_t] id2
        np.ndarray[np.float32_t, ndim=1] np_rg_affs
        np.ndarray[np.uint64_t, ndim=1] np_id1
        np.ndarray[np.uint64_t, ndim=1] np_id2
        size_t i
    
    max_segid = seg.max()
    get_region_graph_average(paff, pseg, max_segid, rg_affs, id1, id2)
    np_rg_affs = np.zeros(rg_affs.size(), np.float32)
    np_id1 = np.zeros(id1.size(), np.uint64)
    np_id2 = np.zeros(id2.size(), np.uint64)
    for 0 <= i < rg_affs.size():
        np_rg_affs[i] = rg_affs[i]
        np_id1[i] = id1[i]
        np_id2[i] = id2[i]
    return (np_rg_affs, np_id1, np_id2)

def zw_merge_segments_with_function(np.ndarray[np.uint64_t, ndim=3] seg,
                                       np.ndarray[np.float32_t, ndim=1] rg_affs,
                                       np.ndarray[np.uint64_t, ndim=1] id1,
                                       np.ndarray[np.uint64_t, ndim=1] id2,
                                       np.ndarray[np.uint64_t, ndim=1] counts,
                                       T_size, T_weight, T_dust, T_merge):
    '''Perform the agglomeration step
    
    :param seg: the segmentation
    :param rg_affs: the affinities from the region graph 3-tuple
                    returned by zw_get_region_graph
    :param id1: the first id from the region graph 3-tuple
    :param id2: the second id from the region graph 3-tuple
    :param counts: the voxel counts returned by zw_find_basins
    :param T_size: the maximum size allowed in the first merge
    :param T_weight: the minimum affinity weight allowed in the first merging
                     step
    :param T_dust: discard objects smaller than this if not merged
    :param T_merge: in the second step, merge if affinities between pairs are
                    greater or equal to this.
    :returns: a two tuple of counts and final region graph. "counts" is a
    one-dimensional numpy array of count per ID. "region graph" is a three
    tuple of numpy arrays: affinity, id1 and id2.
    '''
    cdef:
        PyObject *pseg = <PyObject *>seg
        vector[float] vrg_affs
        vector[uint64_t] vid1
        vector[uint64_t] vid2
        vector[size_t] vcounts
        size_t i
        np.ndarray[np.float32_t, ndim=1] out_rg_affs
        np.ndarray[np.uint64_t, ndim=1] out_id1
        np.ndarray[np.uint64_t, ndim=1] out_id2
        np.ndarray[np.uint64_t, ndim=1] out_counts
        
    assert len(id1) == len(rg_affs)
    assert len(id2) == len(rg_affs)
    for 0 <= i < rg_affs.shape[0]:
        vrg_affs.push_back(rg_affs[i])
        vid1.push_back(id1[i])
        vid2.push_back(id2[i])
    for 0 <= i < counts.shape[0]:
        vcounts.push_back(counts[i])

    merge_segments_with_function(pseg, vrg_affs, vid1, vid2, vcounts,
        T_size, T_weight, T_dust, T_merge);

    out_rg_affs = np.zeros(vrg_affs.size(), np.float32)
    out_id1 = np.zeros(vid1.size(), np.uint64)
    out_id2 = np.zeros(vid2.size(), np.uint64)
    for 0 <= i < vrg_affs.size():
        out_rg_affs[i] = vrg_affs[i]
        out_id1[i] = vid1[i]
        out_id2[i] = vid2[i]
    out_counts = np.zeros(vcounts.size(), np.uint64)
    for 0 <= i < vcounts.size():
        out_counts[i] = vcounts[i]
    return out_counts, (out_rg_affs, out_id1, out_id2)

def zw_mst(np.ndarray[np.float32_t, ndim=1] rg_affs,
           np.ndarray[np.uint64_t, ndim=1] id1,
           np.ndarray[np.uint64_t, ndim=1] id2,
           size_t max_id):
    '''Compute the maximum spanning tree of a region graph
    
    :param rg_affs: region affinities ordered from largest to smallest
    :param id1: the ID of the first region
    :param id2: the ID of the second region
    :returns: a 3-tuple of the reduced rg_affs, id1 and id2
    '''
    cdef:
        vector[float] vrg_affs
        vector[uint64_t] vid1
        vector[uint64_t] vid2
        size_t i
        np.ndarray[np.float32_t, ndim=1] out_rg_affs
        np.ndarray[np.uint64_t, ndim=1] out_id1
        np.ndarray[np.uint64_t, ndim=1] out_id2
    for 0 <= i < rg_affs.shape[0]:
        vrg_affs.push_back(rg_affs[i])
        vid1.push_back(id1[i])
        vid2.push_back(id2[i])
    mst(vrg_affs, vid1, vid2, max_id)
    out_rg_affs = np.zeros(vrg_affs.size(), np.float32)
    out_id1 = np.zeros(vrg_affs.size(), np.uint64)
    out_id2 = np.zeros(vrg_affs.size(), np.uint64)
    for 0 <= i < vrg_affs.size():
        out_rg_affs[i] = vrg_affs[i]
        out_id1[i] = vid1[i]
        out_id2[i] = vid2[i]
    return out_rg_affs, out_id1, out_id2

def zw_do_mapping(np.ndarray[np.uint64_t, ndim=1] id1,
                  np.ndarray[np.uint64_t, ndim=1] id2,
                  np.ndarray[np.uint64_t, ndim=1] counts,
                  size_t max_count):
    '''Find the global mapping of IDs from the region graph
    
    The region graph should be ordered by decreasing affinity and truncated
    at the affinity threshold.
    :param id1: a 1D array of the lefthand side of the two adjacent regions
    :param id2: a 1D array of the righthand side of the two adjacent regions
    :param counts: the number of voxels per ID
    :param max_count: the maximum number of voxels per object
    :returns: a 1D array of the global IDs per local ID
    '''
    cdef:
        vector[uint64_t] vcounts
        vector[uint64_t] vid1
        vector[uint64_t] vid2
        vector[uint64_t] vmapping
        size_t i
        np.ndarray[np.uint64_t, ndim=1] mapping
    
    vmapping.resize(counts.shape[0])
    for 0 <= i < id1.shape[0]:
        vid1.push_back(id1[i])
        vid2.push_back(id2[i])
    for 0 <= i < counts.shape[0]:
        vcounts.push_back(counts[i])
    do_mapping(vid1, vid2, vcounts, vmapping, max_count)
    mapping = np.zeros(vmapping.size(), np.uint64)
    for 0 <= i < vmapping.size():
        mapping[i] = vmapping[i]
    return mapping
    
#-------------- c++ methods --------------------------------------------------------------
cdef extern from "zwatershed.h":
    map[string, list[float]] zw_initial_cpp(size_t dimX, size_t dimY, size_t dimZ, np.float32_t*affs,
                                                np.float32_t thres_low, np.float32_t thres_high)

    map[string, vector[double]] merge_region(size_t dx, size_t dy, size_t dz,
                                               np.float32_t*rgn_graph, int rgn_graph_len, uint64_t*seg,
                                               uint64_t*counts, int counts_len, int thresh, float weight_th, int dust_size, float merge_th)
    void steepest_ascent(PyObject *aff, PyObject *seg, float low, float high)
    void divide_plateaus(PyObject *seg)
    void find_basins(PyObject *seg, vector[uint64_t] &counts)
    void get_region_graph(
        PyObject *aff, PyObject *seg, size_t max_segid, vector[float] &rg_affs,
        vector[uint64_t] &id1, vector[uint64_t] &id2)
    void get_region_graph_average(
        PyObject *aff, PyObject *seg, size_t max_segid, vector[float] &rg_affs,
        vector[uint64_t] &id1, vector[uint64_t] &id2)
    void merge_segments_with_function(
        PyObject *pyseg, vector[float] &rg_affs, vector[uint64_t] &id1,
        vector[uint64_t] &id2, vector[size_t] &counts, size_t size_th,
        float weight_th, size_t lowt, float merge_th)
    void mst(
     vector[float] &rg_affs,
     vector[uint64_t] &id1,
     vector[uint64_t] &id2,
     size_t max_id)
    void do_mapping(
     vector[uint64_t] &id1,
     vector[uint64_t] &id2,
     vector[uint64_t] &counts,
     vector[uint64_t] &mapping,
     uint64_t max_count);
