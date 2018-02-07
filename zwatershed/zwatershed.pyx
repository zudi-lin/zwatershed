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
                  T_dust=600, T_merge=0.5, T_mst=0,
                  seg_save_path=None):
    # aff stats
    affs = np.asfortranarray(np.transpose(affs, (1, 2, 3, 0)))
    dims = affs.shape
    affs_thres = np.percentile(affs, T_aff) if T_aff_relative else T_aff
    print "1. affinity threshold: ", affs_thres

    print "2. get initial seg"
    seg_empty = np.empty((dims[0], dims[1], dims[2]), dtype='uint64')
    map = zwshed_initial(seg_empty, affs, affs_thres[0], affs_thres[1])

    # get initial rg
    cdef np.ndarray[uint64_t, ndim=1] seg_in = map['seg']
    cdef np.ndarray[uint64_t, ndim=1] counts_out = map['counts']
    cdef np.ndarray[np.float32_t, ndim=2] rgn_graph = map['rg']

    # get segs, stats
    T_threshes.sort()
    segs = [None]*len(T_threshes)
    for i in range(len(T_threshes)):
        print "3. do thres: ", T_threshes[i], T_dust
        if(len(rgn_graph) > 0):
            map = merge_region(
                dims[0], dims[1], dims[2], &rgn_graph[0, 0],
                rgn_graph.shape[0], &seg_in[0], &counts_out[0], 
                len(map['counts']), T_threshes[i], affs_thres[2], T_dust, T_merge, T_mst)
        seg = np.array(map['seg'], dtype='uint64').reshape((dims[2], dims[1], dims[0])).transpose(2, 1, 0)
        graph = np.array(map['rg'], dtype='float32')
        counts_out = np.array(map['counts'], dtype='uint64')
        counts_len = len(counts_out)
        seg_in = np.array(map['seg'], dtype='uint64')
        rgn_graph = graph.reshape(len(graph) / 3, 3)
        return rgn_graph

        if seg_save_path is None:
            segs[i] = seg.copy()
        else:
            f = h5py.File(seg_save_path + '_' + str(T_threshes[i]) + '.h5', 'w')
            ds = f.create_dataset('stack', seg.shape, compression="gzip", dtype=np.uint64)
            ds[:] = seg.astype(np.uint32)
            f.close()
        print "\t number of regions: "+ str(len(np.unique(seg)))
    if seg_save_path is None:
        return segs

#################################################
# auxilary function for debug purpose

def zwshed_initial(np.ndarray[uint64_t, ndim=3] seg, np.ndarray[np.float32_t, ndim=4] affs, affs_low, affs_high):
    cdef np.ndarray[uint64_t, ndim=1] counts = np.empty(1, dtype='uint64')
    dims = affs.shape
    map = zwshed_initial_c(dims[0], dims[1], dims[2], &affs[0, 0, 0, 0], float(affs_low), float(affs_high))
    graph = np.array(map['rg'], dtype='float32')
    return {'rg': graph.reshape(len(graph) / 3, 3), 'seg': np.array(map['seg'], dtype='uint64'),
            'counts': np.array(map['counts'], dtype='uint64')}

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
                        np.ndarray[np.uint64_t, ndim=3] seg,
                        max_segid):
    '''Return the initial region graph using max edge affinity.
    
    :param aff: the affinity predictions - an array of x, y, z, c where c == 0
                is the affinity prediction for x, c == 1 is the affinity
                prediction for y and c == 2 is the affinity prediction for z
    :param seg: the segmentation after finding basins
    :param max_segid: the maximum ID in seg
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
                                np.ndarray[np.uint64_t, ndim=3] seg,
                                max_segid):
    '''Return the initial region graph using average edge affinity
    
    :param aff: the affinity predictions - an array of x, y, z, c where c == 0
                is the affinity prediction for x, c == 1 is the affinity
                prediction for y and c == 2 is the affinity prediction for z
    :param seg: the segmentation after finding basins
    :param max_segid: the maximum ID in seg
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
    
    get_region_graph_average(paff, pseg, max_segid, rg_affs, id1, id2)
    np_rg_affs = np.zeros(rg_affs.size(), np.float32)
    np_id1 = np.zeros(id1.size(), np.uint64)
    np_id2 = np.zeros(id2.size(), np.uint64)
    for 0 <= i < rg_affs.size():
        np_rg_affs[i] = rg_affs[i]
        np_id1[i] = id1[i]
        np_id2[i] = id2[i]
    return (np_rg_affs, np_id1, np_id2)

def zw_merge_segments_with_function_dw(np.ndarray[np.uint64_t, ndim=3] seg,
                                       np.ndarray[np.float32_t, ndim=1] rg_affs,
                                       np.ndarray[np.uint64_t, ndim=1] id1,
                                       np.ndarray[np.uint64_t, ndim=1] id2,
                                       np.ndarray[np.uint64_t, ndim=1] counts,
                                       T_size, T_weight, T_dust, T_merge, T_mst):
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
    :param T_mst: if do MST
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

    merge_segments_with_function_dw(pseg, vrg_affs, vid1, vid2, vcounts,
        T_size, T_weight, T_dust, T_merge, T_mst);

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
    
    
#-------------- c++ methods --------------------------------------------------------------
cdef extern from "zwatershed.h":
    map[string, list[float]] zwshed_initial_c(size_t dimX, size_t dimY, size_t dimZ, np.float32_t*affs,
                                                np.float32_t thres_low, np.float32_t thres_high)

    map[string, vector[double]] merge_region(size_t dx, size_t dy, size_t dz,
                                               np.float32_t*rgn_graph, int rgn_graph_len, uint64_t*seg,
                                               uint64_t*counts, int counts_len, int thresh, float weight_th, int dust_size, float merge_th, int do_mst)
    void steepest_ascent(PyObject *aff, PyObject *seg, float low, float high)
    void divide_plateaus(PyObject *seg)
    void find_basins(PyObject *seg, vector[uint64_t] &counts)
    void get_region_graph(
        PyObject *aff, PyObject *seg, size_t max_segid, vector[float] &rg_affs,
        vector[uint64_t] &id1, vector[uint64_t] &id2)
    void get_region_graph_average(
        PyObject *aff, PyObject *seg, size_t max_segid, vector[float] &rg_affs,
        vector[uint64_t] &id1, vector[uint64_t] &id2)
    void merge_segments_with_function_dw(
        PyObject *pyseg, vector[float] &rg_affs, vector[uint64_t] &id1,
        vector[uint64_t] &id2, vector[size_t] &counts, size_t size_th,
        float weight_th, size_t lowt, float merge_th, int do_mst)


