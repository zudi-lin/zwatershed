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
def zwatershed_dw(np.ndarray[np.float32_t, ndim=4] affs, 
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
    map = zwshed_initial_dw(seg_empty, affs, affs_thres[0], affs_thres[1])

    # get initial rg
    cdef np.ndarray[uint64_t, ndim=1] seg_in = map['seg']
    cdef np.ndarray[uint64_t, ndim=1] counts_out = map['counts']
    cdef np.ndarray[np.float32_t, ndim=2] rgn_graph = map['rg']

    # get segs, stats
    threshes.sort()
    for i in range(len(threshes)):
        print "3. do thres: ", threshes[i], T_dust
        if(len(rgn_graph) > 0):
            map = merge_no_stats_dw(
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

def zwshed_initial_dw(np.ndarray[uint64_t, ndim=3] seg, np.ndarray[np.float32_t, ndim=4] affs, affs_low, affs_high):
    cdef np.ndarray[uint64_t, ndim=1] counts = np.empty(1, dtype='uint64')
    dims = affs.shape
    map = zwshed_initial_c_dw(dims[0], dims[1], dims[2], &affs[0, 0, 0, 0], float(affs_low), float(affs_high))
    graph = np.array(map['rg'], dtype='float32')
    return {'rg': graph.reshape(len(graph) / 3, 3), 'seg': np.array(map['seg'], dtype='uint64'),
            'counts': np.array(map['counts'], dtype='uint64')}


def zwatershed_and_metrics(gt, affs, threshes, save_threshes):
    threshes.sort()
    return zwshed_with_stats(gt, affs, threshes, save_threshes, h5=0)

def zwatershed_and_metrics_h5(gt, affs, threshes, save_threshes, seg_save_path):
    threshes.sort()
    return zwshed_with_stats(gt, affs, threshes, save_threshes, h5=1, seg_save_path=seg_save_path)


def zwatershed(affs, threshes):
    threshes.sort()
    return zwshed_no_stats(affs, threshes, threshes, h5=0)

def zwatershed_h5(affs, threshes, seg_save_path):
    threshes.sort()
    zwshed_no_stats(affs, threshes, threshes, h5=1, seg_save_path=seg_save_path)

def zwatershed_and_metrics_arb(gt, node1, node2, edgeWeight, threshes, save_threshes):
    threshes.sort()
    return zwshed_with_stats_arb(gt, node1, node2, edgeWeight, threshes, save_threshes, h5=0)

def zwatershed_and_metrics_h5_arb(gt, node1, node2, edgeWeight, threshes, save_threshes, seg_save_path):
    threshes.sort()
    return zwshed_with_stats_arb(gt, node1, node2, edgeWeight, threshes, save_threshes, h5=1,
                                 seg_save_path=seg_save_path)

def zwatershed_arb(seg_shape, node1, node2, edgeWeight, save_threshes):
    save_threshes.sort()
    return zwshed_no_stats_arb(seg_shape, node1, node2, edgeWeight, save_threshes, save_threshes, h5=0)

def zwatershed_h5_arb(seg_shape, node1, node2, edgeWeight, save_threshes, seg_save_path):
    save_threshes.sort()
    return zwshed_no_stats_arb(seg_shape, node1, node2, edgeWeight, save_threshes, save_threshes, h5=1,
                               seg_save_path=seg_save_path)
                               
def zwatershed_basic_h5(affs, seg_save_path):
    zwatershed_basic_h5(affs, seg_save_path=seg_save_path)

#-------------- helper methods --------------------------------------------------------------
def zwatershed_basic_h5(np.ndarray[np.float32_t, ndim=4] affs, seg_save_path="NULL/"):
    makedirs(seg_save_path)

    # get initial seg,rg
    affs = np.asfortranarray(np.transpose(affs, (1, 2, 3, 0)))
    dims = affs.shape
    seg_empty = np.empty((dims[0], dims[1], dims[2]), dtype='uint64')
    map = zwshed_initial(seg_empty, affs)
    counts = map['counts']
    rg = map['rg']
    f = h5py.File(seg_save_path + 'basic.h5', 'w')
    f["seg"] = seg = np.array(map['seg'], dtype='uint64').reshape((dims[2], dims[1], dims[0])).transpose(2, 1, 0)
    f["counts"]=counts
    f["rg"]=rg
    f.close()

def zwshed_with_stats(np.ndarray[uint64_t, ndim=3] gt, np.ndarray[np.float32_t, ndim=4] affs, threshes, save_threshes,
                      int h5, seg_save_path="NULL/"):
    if h5:
        makedirs(seg_save_path)

    # get initial seg,rg
    affs = np.asfortranarray(np.transpose(affs, (1, 2, 3, 0)))
    gt = np.array(gt, order='F')
    map = zwshed_initial(gt, affs)
    cdef np.ndarray[uint64_t, ndim=1] seg_in = map['seg']
    cdef np.ndarray[uint64_t, ndim=1] counts_out = map['counts']
    cdef np.ndarray[np.float32_t, ndim=2] rgn_graph = map['rg']

    counts_len = len(map['counts'])
    dims = affs.shape

    # get segs, stats
    segs, splits, merges, info_splits, info_merges = [], [], [], [], []
    for i in range(len(threshes)):
        if(len(rgn_graph) > 0):
            map = merge_with_stats(dims[0], dims[1], dims[2], &gt[0, 0, 0], &rgn_graph[0, 0],
                               rgn_graph.shape[0], &seg_in[0], &counts_out[0], counts_len, threshes[i])
        seg = np.array(map['seg'], dtype='uint64').reshape((dims[2], dims[1], dims[0])).transpose(2, 1, 0)
        graph = np.array(map['rg'], dtype='float32')
        counts_out = np.array(map['counts'], dtype='uint64')
        counts_len = len(counts_out)
        seg_in = np.array(map['seg'], dtype='uint64')
        rgn_graph = graph.reshape(len(graph) / 3, 3)
        if threshes[i] in save_threshes:
            if h5:
                f = h5py.File(seg_save_path + 'seg_' + str(threshes[i]) + '.h5', 'w')
                f["main"] = seg
                f.close()
            else:
                segs.append(seg)
        splits = splits + [map['stats'][0]]
        merges = merges + [map['stats'][1]]
        info_splits = info_splits + [map['stats'][2]]
        info_merges = info_merges + [map['stats'][3]]
    max_f_score = 2 / (1 / splits[0] + 1 / merges[0])
    max_v_info = 2 / (1 / info_splits[0] + 1 / info_merges[0])
    for j in range(len(splits)):
        f_score = 2 / (1 / splits[j] + 1 / merges[j])
        if f_score > max_f_score:
            max_f_score = f_score
        info_score = 2 / (1 / info_splits[j] + 1 / info_merges[j])
        if info_score > max_v_info:
            max_v_info = info_score
        
    returnMap = {'V_Rand': max_f_score, 'V_Rand_split': splits, 'V_Rand_merge': merges, 'V_Info_split':info_splits, 'V_Info_merge':info_merges, 'V_Info':max_v_info}
    if h5:
        return returnMap
    else:
        return segs, returnMap


def zwshed_no_stats(np.ndarray[np.float32_t, ndim=4] affs, threshes, save_threshes, int h5,
                    seg_save_path="NULL/"):
    if h5:
        makedirs(seg_save_path)

    # get initial seg,rg
    affs = np.asfortranarray(np.transpose(affs, (1, 2, 3, 0)))
    dims = affs.shape
    seg_empty = np.empty((dims[0], dims[1], dims[2]), dtype='uint64')
    map = zwshed_initial(seg_empty, affs)
    cdef np.ndarray[uint64_t, ndim=1] seg_in = map['seg']
    cdef np.ndarray[uint64_t, ndim=1] counts_out = map['counts']
    cdef np.ndarray[np.float32_t, ndim=2] rgn_graph = map['rg']
    segs = []

    # get segs, stats
    for i in range(len(threshes)):
        if(len(rgn_graph) > 0):
            map = merge_no_stats(dims[0], dims[1], dims[2], &rgn_graph[0, 0],
                             rgn_graph.shape[0], &seg_in[0], &counts_out[0], len(map['counts']), threshes[i])
        seg = np.array(map['seg'], dtype='uint64').reshape((dims[2], dims[1], dims[0])).transpose(2, 1, 0)
        graph = np.array(map['rg'], dtype='float32')
        counts_out = np.array(map['counts'], dtype='uint64')
        counts_len = len(counts_out)
        seg_in = np.array(map['seg'], dtype='uint64')
        rgn_graph = graph.reshape(len(graph) / 3, 3)
        if threshes[i] in save_threshes:
            if h5:
                f = h5py.File(seg_save_path + 'seg_' + str(threshes[i]) + '.h5', 'w')
                f["main"] = seg
                f.close()
            else:
                segs.append(seg)
    if not h5:
        return segs

def zwshed_initial(np.ndarray[uint64_t, ndim=3] seg, np.ndarray[np.float32_t, ndim=4] affs):
    cdef np.ndarray[uint64_t, ndim=1] counts = np.empty(1, dtype='uint64')
    dims = affs.shape
    map = zwshed_initial_c(dims[0], dims[1], dims[2], &affs[0, 0, 0, 0])
    graph = np.array(map['rg'], dtype='float32')
    return {'rg': graph.reshape(len(graph) / 3, 3), 'seg': np.array(map['seg'], dtype='uint64'),
            'counts': np.array(map['counts'], dtype='uint64')}

def makedirs(seg_save_path):
    if not seg_save_path.endswith("/"):
        seg_save_path += "/"
    if not os.path.exists(seg_save_path):
        os.makedirs(seg_save_path)

# arb methods ---------------
def zwshed_with_stats_arb(np.ndarray[uint64_t, ndim=3] gt, np.ndarray[uint64_t, ndim=1] node1,
                          np.ndarray[uint64_t, ndim=1] node2, np.ndarray[float, ndim=1] edgeWeight, threshes,
                          save_threshes,
                          int h5, seg_save_path="NULL/"):
    if h5:
        makedirs(seg_save_path)

    # get initial seg,rg
    gt = np.array(gt, order='C')
    n_edge = node1.size
    map = zwshed_initial_arb(gt, node1, node2, n_edge, edgeWeight)
    cdef np.ndarray[uint64_t, ndim=1] seg_in = map['seg']
    cdef np.ndarray[uint64_t, ndim=1] counts_out = map['counts']
    cdef np.ndarray[np.float32_t, ndim=2] rgn_graph = map['rg']

    counts_len = len(map['counts'])
    dims = gt.shape

    # get segs, stats
    segs, splits, merges = [], [], []
    for i in range(len(threshes)):
        if(len(rgn_graph) > 0):
            map = merge_with_stats_arb(dims[0], dims[1], dims[2], &gt[0, 0, 0], &rgn_graph[0, 0],
                                   rgn_graph.shape[0], &seg_in[0], &counts_out[0], counts_len, threshes[i])
        seg = np.array(map['seg'], dtype='uint64').reshape((dims[0], dims[1], dims[2]))
        graph = np.array(map['rg'], dtype='float32')
        counts_out = np.array(map['counts'], dtype='uint64')
        counts_len = len(counts_out)
        seg_in = np.array(map['seg'], dtype='uint64')
        rgn_graph = graph.reshape(len(graph) / 3, 3)
        if threshes[i] in save_threshes:
            if h5:
                f = h5py.File(seg_save_path + 'seg_' + str(threshes[i]) + '.h5', 'w')
                f["main"] = seg
                f.close()
            else:
                segs.append(seg)
        splits = splits + [map['stats'][0]]
        merges = merges + [map['stats'][1]]
    max_f_score = 2 / (1 / splits[0] + 1 / merges[0])
    for j in range(len(splits)):
        f_score = 2 / (1 / splits[j] + 1 / merges[j])
        if f_score > max_f_score:
            max_f_score = f_score
    returnMap = {'V_Rand': max_f_score, 'V_Rand_split': splits, 'V_Rand_merge': merges}
    if h5:
        return returnMap
    else:
        return segs, returnMap

def zwshed_no_stats_arb(dims, np.ndarray[uint64_t, ndim=1] node1,
                        np.ndarray[uint64_t, ndim=1] node2, np.ndarray[float, ndim=1] edgeWeight, threshes,
                        save_threshes,
                        int h5, seg_save_path="NULL/"):
    if h5:
        makedirs(seg_save_path)

    # get initial seg,rg
    n_edge = node1.size
    seg_empty = np.zeros(dims,dtype='uint64')
    map = zwshed_initial_arb(seg_empty, node1, node2, n_edge, edgeWeight)
    cdef np.ndarray[uint64_t, ndim=1] seg_in = map['seg']
    cdef np.ndarray[uint64_t, ndim=1] counts_out = map['counts']
    cdef np.ndarray[np.float32_t, ndim=2] rgn_graph = map['rg']
    counts_len = len(map['counts'])

    # get segs, stats
    segs = []
    for i in range(len(threshes)):
        if(len(rgn_graph) > 0):
            map = merge_no_stats_arb(dims[0], dims[1], dims[2], &rgn_graph[0, 0],
                                 rgn_graph.shape[0], &seg_in[0], &counts_out[0], counts_len, threshes[i])
        seg = np.array(map['seg'], dtype='uint64').reshape((dims[0], dims[1], dims[2]))
        graph = np.array(map['rg'], dtype='float32')
        counts_out = np.array(map['counts'], dtype='uint64')
        counts_len = len(counts_out)
        seg_in = np.array(map['seg'], dtype='uint64')
        rgn_graph = graph.reshape(len(graph) / 3, 3)
        if threshes[i] in save_threshes:
            if h5:
                f = h5py.File(seg_save_path + 'seg_' + str(threshes[i]) + '.h5', 'w')
                f["main"] = seg
                f.close()
            else:
                segs.append(seg)
    if not h5:
        return segs

def zwshed_initial_arb(np.ndarray[uint64_t, ndim=3] seg, np.ndarray[uint64_t, ndim=1] node1,
                       np.ndarray[uint64_t, ndim=1] node2, int n_edge, np.ndarray[float, ndim=1] edgeWeight):
    cdef np.ndarray[uint64_t, ndim=1] counts = np.empty(1, dtype='uint64')
    dims = seg.shape
    map = zwshed_initial_c_arb(dims[0], dims[1], dims[2], &node1[0], &node2[0], &edgeWeight[0], n_edge)
    graph = np.array(map['rg'], dtype='float32')
    return {'rg': graph.reshape(len(graph) / 3, 3), 'seg': np.array(map['seg'], dtype='uint64'),
            'counts': np.array(map['counts'], dtype='uint64')}



#-------------- c++ methods --------------------------------------------------------------
cdef extern from "zwatershed.h":
    map[string, list[float]] zwshed_initial_c_dw(size_t dimX, size_t dimY, size_t dimZ, np.float32_t*affs,
                                                np.float32_t thres_low, np.float32_t thres_high)

    map[string, vector[double]] merge_no_stats_dw(size_t dx, size_t dy, size_t dz,
                                               np.float32_t*rgn_graph, int rgn_graph_len, uint64_t*seg,
                                               uint64_t*counts, int counts_len, int thresh, float weight_th, int dust_size, float merge_th)
    # ------------
    map[string, list[float]] zwshed_initial_c(size_t dimX, size_t dimY, size_t dimZ, np.float32_t*affs)

    map[string, vector[double]] merge_with_stats(size_t dx, size_t dy, size_t dz, np.uint64_t*gt,
                                                 np.float32_t*rgn_graph, int rgn_graph_len, uint64_t*seg,
                                                 uint64_t*counts, int counts_len, int thresh)
    map[string, vector[double]] merge_no_stats(size_t dx, size_t dy, size_t dz,
                                               np.float32_t*rgn_graph, int rgn_graph_len, uint64_t*seg,
                                               uint64_t*counts, int counts_len, int thresh)
    map[string, list[float]] zwshed_initial_c_arb(size_t dimX, size_t dimY, size_t dimZ, uint64_t*node1,
                                                  uint64_t*node2, float*edgeWeight, int n_edge)
    map[string, vector[double]] merge_with_stats_arb(size_t dx, size_t dy, size_t dz, np.uint64_t*gt,
                                                     np.float32_t*rgn_graph, int rgn_graph_len, uint64_t*seg,
                                                     uint64_t*counts, int counts_len, int thresh)
    map[string, vector[double]] merge_no_stats_arb(size_t dx, size_t dy, size_t dz,
                                                   np.float32_t*rgn_graph, int rgn_graph_len, uint64_t*seg,
                                                   uint64_t*counts, int counts_len, int thresh)

    void steepest_ascent(PyObject *aff, PyObject *seg, float low, float high)
    void divide_plateaus(PyObject *seg)
    void find_basins(PyObject *seg, vector[uint64_t] &counts)

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