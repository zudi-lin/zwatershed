#ifndef TESTLIB_H
#define TESTLIB_H

#include <iostream>
#include <list>
#include <string>
#include <map>
#include <vector>
#include <utility>


std::map<std::string,std::list<float>> zwshed_initial_c_dw(const size_t dx, const size_t dy, const size_t dz, float* affs,
                                                float thres_low, float thres_high);
std::map<std::string,std::vector<double>> merge_no_stats_dw(size_t dimX, size_t dimY, size_t dimZ, float * rgn_graph,
                                        int rgn_graph_len, uint64_t * seg_in, uint64_t*counts, int counts_len, 
                                        int thresh, float T_aff_merge, int T_dust, float T_merge);
// ------------

std::map<std::string,std::list<float>> zwshed_initial_c(const size_t dx, const size_t dy, const size_t dz, float* affs);
std::map<std::string,std::vector<double>> merge_with_stats(size_t dx,size_t dy, size_t dz, uint64_t * gt, float * rgn_graph,
                                        int rgn_graph_len, uint64_t * seg_in, uint64_t*counts, int counts_len, int thresh);

std::map<std::string,std::vector<double>> merge_no_stats(size_t dimX, size_t dimY, size_t dimZ, float * rgn_graph,
                                        int rgn_graph_len, uint64_t * seg_in, uint64_t*counts, int counts_len, int thresh);

std::map<std::string,std::list<float>> zwshed_initial_c_arb(const size_t dx, const size_t dy, const size_t dz, const uint64_t*node1,
                                               const uint64_t*node2, const float*edgeWeight, const int n_edge);

std::map<std::string,std::vector<double>> merge_with_stats_arb(size_t dx,size_t dy, size_t dz, uint64_t * gt, float * rgn_graph,
                                        int rgn_graph_len, uint64_t * seg_in, uint64_t*counts, int counts_len, int thresh);

std::map<std::string,std::vector<double>> merge_no_stats_arb(size_t dx,size_t dy, size_t dz, float * rgn_graph,
                                        int rgn_graph_len, uint64_t * seg_in, uint64_t*counts, int counts_len, int thresh);

#endif
