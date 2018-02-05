/* Connected components
 * developed and maintained by Srinivas C. Turaga <sturaga@mit.edu>
 * do not distribute without permission.
 */
#include "zwatershed.h"
//#pragma once
#include "zwatershed_util/agglomeration.hpp"
#include "zwatershed_util/region_graph.hpp"
#include "zwatershed_util/basic_watershed.hpp"
#include "zwatershed_util/limit_functions.hpp"
#include "zwatershed_util/types.hpp"
#include "zwatershed_util/main_helper.hpp"


#include <memory>
#include <type_traits>

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstddef>
#include <cstdint>
#include <queue>
#include <vector>
#include <algorithm>
#include <tuple>
#include <map>
#include <list>
#include <set>
#include <vector>
#include <chrono>
#include <fstream>
#include <string>
#include <boost/make_shared.hpp>
using namespace std;
// these values based on 5% at iter = 10000
double LOW=  0.0001;
double HIGH= 0.9999;
bool RECREATE_RG = true;

std::map<std::string,std::list<float>> zwshed_initial_c(const size_t dimX, const size_t dimY, const size_t dimZ, float* affs,
                                                           float thres_low, float thres_high)
{
    std::cout << "calculating basic watershed..." << std::endl;

    // read data
    volume_ptr<uint64_t> seg_ref;
    std::vector<std::size_t> counts_ref;
    affinity_graph_ptr<float> aff(new affinity_graph<float>
                              (boost::extents[dimX][dimY][dimZ][3],
                               boost::fortran_storage_order()));
    for(size_t i=0;i<dimX*dimY*dimZ*3;i++)
        aff->data()[i] = affs[i];
    std::tie(seg_ref , counts_ref) = watershed<uint64_t>(aff, thres_low, thres_high);


    // calculate region graph
    std::cout << "calculating rgn graph..." << std::endl;
    auto rg = get_region_graph(*aff, *seg_ref , counts_ref.size()-1);
    std::cout << "Finished calculating region graph" << std::endl;
    // save and return
    std::map<std::string,std::list<float>> returnMap;
    std::list<float> empty;
    returnMap["rg"] = empty;
    std::list<float> &rg_data = returnMap["rg"];
    for ( const auto& e: *rg ){
        rg_data.push_back(std::get<1>(e));
        rg_data.push_back(std::get<2>(e));
        rg_data.push_back(std::get<0>(e));
    }
    std::cout << "Copied region graph" << std::endl;
    returnMap["seg"] = empty;
    std::list<float> &seg_data = returnMap["seg"];
    returnMap["counts"] = empty;
    std::list<float> &counts_data = returnMap["counts"];
    for(size_t i=0;i<dimX*dimY*dimZ;i++)
        seg_data.push_back(seg_ref->data()[i]);
    std::cout << "copied segmentation" << std::endl;
    for (const auto& x:counts_ref)
        counts_data.push_back(x);
    std::cout << "copied counts" << std::endl;
    std::cout << "Returning from zwshed_initial_c_dw" << std::endl;
    return returnMap;
 }
std::map<std::string,std::vector<double>> merge_no_stats(
    size_t dimX, size_t dimY, size_t dimZ, float * rgn_graph, int rgn_graph_len,
    uint64_t * seg_in, uint64_t*counts_in, int counts_len, int thresh, 
    float T_aff_merge, int T_dust, float T_merge) {
    std::cout << "evaluating..." << std::endl;

    // read data
    volume_ptr<uint64_t> seg(new volume<uint64_t> (boost::extents[dimX][dimY][dimZ], boost::fortran_storage_order()));
    std::vector<std::size_t> counts = * new std::vector<std::size_t>();
    region_graph_ptr<uint64_t,float> rg( new region_graph<uint64_t,float> );
    for(size_t i=0;i<dimX*dimY*dimZ;i++)
        seg->data()[i] = seg_in[i];
    for(int i=0;i<counts_len;i++)
        counts.push_back(counts_in[i]);
    for(int i=0;i<rgn_graph_len;i++)
        (*rg).emplace_back(rgn_graph[i*3+2],rgn_graph[i*3],rgn_graph[i*3+1]);

    // merge
    std::cout << "thresh: " << thresh << "\n";
    rg = merge_segments_with_function_dw(
	  *seg, *rg, counts, thresh, T_aff_merge, T_dust , T_merge);

	// save and return
	std::map<std::string,std::vector<double>> returnMap;
    std::vector<double> seg_vector;
    std::vector<double> rg_data; // = * (new std::list<float>());
    std::vector<double> counts_data; // = * (new std::list<float>());
    for(size_t i=0;i<dimX*dimY*dimZ;i++)
        seg_vector.push_back(((double)(seg->data()[i])));
    for ( const auto& e: *rg ){
        rg_data.push_back(std::get<1>(e));
        rg_data.push_back(std::get<2>(e));
        rg_data.push_back(std::get<0>(e));
    }
    for (const auto& x:counts)
        counts_data.push_back(x);
    returnMap["seg"] = seg_vector;
    returnMap["rg"]=rg_data;
    returnMap["counts"] = counts_data;
    return returnMap;
 }

/*
 * Affinity is organized into a 4-d array. The first index is the
 * Z / Y / X affinity in question. This index is of length = 3.
 * The remaining indices are the Z, Y and X coordinates.
 */
class affinity_t: public boost::multi_array_ref<float, 4> {
    public:
	typedef boost::multi_array_ref<float, 4> super;
	affinity_t(PyArrayObject *a):super(
	    (float *)PyArray_DATA(a), 
	    boost::extents[PyArray_DIMS(a)[0]]
			  [PyArray_DIMS(a)[1]]
			  [PyArray_DIMS(a)[2]]
			  [PyArray_DIMS(a)[3]]) {
	    for (size_t i=0; i < 4; i++) {
		stride_list_[i] = PyArray_STRIDE(a, i) / sizeof(float);
	    }
	}
};
/*
 * The seeds and segmentation are hardcoded as uint32 for use by Python.
 */
class segmentation_t: public boost::multi_array_ref<uint64_t, 3> {
    public:
	typedef boost::multi_array_ref<uint64_t, 3> super;
	segmentation_t(PyArrayObject *a):super(
	    (uint64_t *)PyArray_DATA(a), 
	    boost::extents[PyArray_DIMS(a)[0]]
			 [PyArray_DIMS(a)[1]]
			 [PyArray_DIMS(a)[2]]) {
	    for (size_t i=0; i < 3; i++) {
		stride_list_[i] = PyArray_STRIDE(a, i) / sizeof(uint64_t);
	    }
	}
};

void steepest_ascent(PyObject *pyaff, PyObject *pyseg, float low, float high) {
    affinity_t aff((PyArrayObject *)pyaff);
    segmentation_t seg((PyArrayObject *)pyseg);
    steepestascent(aff, seg, low, high);
}

void divide_plateaus(PyObject *pyseg) {
    segmentation_t seg((PyArrayObject *)pyseg);
    divideplateaus(seg);
}

void find_basins(PyObject *pyseg, std::vector<uint64_t> &counts) {
    segmentation_t seg((PyArrayObject *)pyseg);
    findbasins(seg, counts);
}

void rg_to_vectors(const region_graph<uint64_t, float> &rg,
		   std::vector<float> &rg_affs, 
		   std::vector<uint64_t> &id1, std::vector<uint64_t> &id2) {
    rg_affs.clear();
    id1.clear();
    id2.clear();
    for (auto e:rg) {
	rg_affs.push_back(std::get<0>(e));
	id1.push_back(std::get<1>(e));
	id2.push_back(std::get<2>(e));
    }
}

void rg_from_vectors(region_graph<uint64_t, float> &rg,
		     const std::vector<float> &rg_affs, 
		     const std::vector<uint64_t> &id1, 
		     const std::vector<uint64_t> &id2) {
    size_t i;
    rg.clear();
    for (i=0; i<rg_affs.size(); i++) {
	rg.emplace_back(rg_affs[i], id1[i], id2[i]);
    }
}

void get_region_graph(
    PyObject *pyaff, PyObject *pyseg, std::size_t max_segid,
    std::vector<float> &rg_affs, 
    std::vector<uint64_t> &id1, std::vector<uint64_t> &id2) {
	
    affinity_t aff((PyArrayObject *)pyaff);
    segmentation_t seg((PyArrayObject *)pyseg);
    auto rg_ptr = get_region_graph(aff, seg, max_segid);
    rg_to_vectors(*rg_ptr, rg_affs, id1, id2);
}

void get_region_graph_average(
    PyObject *pyaff, PyObject *pyseg, std::size_t max_segid,
    std::vector<float> &rg_affs, 
    std::vector<uint64_t> &id1, std::vector<uint64_t> &id2) {
	
    affinity_t aff((PyArrayObject *)pyaff);
    segmentation_t seg((PyArrayObject *)pyseg);
    auto rg_ptr = get_region_graph_average(aff, seg, max_segid);
    rg_to_vectors(*rg_ptr, rg_affs, id1, id2);
}

void merge_segments_with_function_dw(
     PyObject *pyseg,
     std::vector<float> &rg_affs,
     std::vector<uint64_t> &id1,
     std::vector<uint64_t> &id2,
     std::vector<std::size_t> &counts,
     const size_t size_th,
     const float weight_th,
     const size_t lowt,
     const float merge_th) {
    segmentation_t seg((PyArrayObject *)pyseg);
    region_graph<uint64_t, float> rg;
    
    rg_from_vectors(rg, rg_affs, id1, id2);
    auto rg_ptr = merge_segments_with_function_dw(
	seg, rg, counts, size_th, weight_th, lowt, merge_th);
    rg_to_vectors(*rg_ptr, rg_affs, id1, id2);
}
     
