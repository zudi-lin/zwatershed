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
                  threshes, T_aff=[0.01,0.8,0.2], 
                  T_dust=600, 
                  T_merge=0.5,
                  seg_save_path='./', T_aff_relative=True):
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
    threshes.sort()
    for i in range(len(threshes)):
        print "3. do thres: ", threshes[i], T_dust
        if(len(rgn_graph) > 0):
            map = merge_no_stats(
                dims[0], dims[1], dims[2], &rgn_graph[0, 0],
                rgn_graph.shape[0], &seg_in[0], &counts_out[0], 
                len(map['counts']), threshes[i], affs_thres[2], T_dust, T_merge)
        seg = np.array(map['seg'], dtype='uint64').reshape((dims[2], dims[1], dims[0])).transpose(2, 1, 0)
        graph = np.array(map['rg'], dtype='float32')
        counts_out = np.array(map['counts'], dtype='uint64')
        counts_len = len(counts_out)
        seg_in = np.array(map['seg'], dtype='uint64')
        rgn_graph = graph.reshape(len(graph) / 3, 3)
        f = h5py.File(seg_save_path + '_' + str(threshes[i]) + '.h5', 'w')
        ds = f.create_dataset('stack', seg.shape, compression="gzip", dtype=np.uint32)
        ds[:] = seg.astype(np.uint32)
        f.close()
        print "\t number of regions: "+ str(len(np.unique(seg)))

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

#-------------- c++ methods --------------------------------------------------------------
cdef extern from "zwatershed.h":
    map[string, list[float]] zwshed_initial_c(size_t dimX, size_t dimY, size_t dimZ, np.float32_t*affs,
                                                np.float32_t thres_low, np.float32_t thres_high)

    map[string, vector[double]] merge_no_stats(size_t dx, size_t dy, size_t dz,
                                               np.float32_t*rgn_graph, int rgn_graph_len, uint64_t*seg,
                                               uint64_t*counts, int counts_len, int thresh, float weight_th, int dust_size, float merge_th)
    void steepest_ascent(PyObject *aff, PyObject *seg, float low, float high)
    void divide_plateaus(PyObject *seg)
    void find_basins(PyObject *seg, vector[uint64_t] &counts)


