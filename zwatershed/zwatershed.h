#ifndef TESTLIB_H
#define TESTLIB_H

#include <iostream>
#include <list>
#include <string>
#include <map>
#include <tuple>
#include <vector>
#include <utility>
#include <Python.h>
#include "numpy/arrayobject.h"

std::map<std::string,std::list<float>> zw_initial_cpp(const size_t dx, const size_t dy, const size_t dz, float* affs,
                                                float thres_low, float thres_high);
std::map<std::string,std::vector<double>> merge_region(size_t dimX, size_t dimY, size_t dimZ, float * rgn_graph,
                                        int rgn_graph_len, uint64_t * seg_in, uint64_t*counts, int counts_len, 
                                        int thresh, float T_aff_merge, int T_dust, float T_merge);

void steepest_ascent(PyObject *aff, PyObject *seg, float low, float high);
void divide_plateaus(PyObject *seg);
void find_basins(PyObject *seg, std::vector<uint64_t> &counts);
void get_region_graph(
    PyObject *aff, PyObject *seg, size_t max_segid,
    std::vector<float> &rg_affs, 
    std::vector<uint64_t> &id1, std::vector<uint64_t> &id2);
void get_region_graph_average(
    PyObject *aff, PyObject *seg, size_t max_segid,
    std::vector<float> &rg_affs, 
    std::vector<uint64_t> &id1, std::vector<uint64_t> &id2);
void merge_segments_with_function(
     PyObject *pyseg,
     std::vector<float> &rg_affs,
     std::vector<uint64_t> &id1,
     std::vector<uint64_t> &id2,
     std::vector<std::size_t> &counts,
     const size_t size_th,
     const float weight_th,
     const size_t lowt,
     const float merge_th);
#endif
