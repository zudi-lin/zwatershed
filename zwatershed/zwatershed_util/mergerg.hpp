/* mergerg - merge the maximum spanning tree region graph
 *
 * Adapted from mergrg! in https://github.com/seung-lab/segment.jl
 */
#pragma once
#include "types.hpp"
#include <map>
#include <tuple>
#include <set>
#include <iostream>
/*
 * mergerg - merge all nodes in the region graph above a threshold
 *
 * Template parameters:
 * ID - the segmentation ID type
 * F - the type of the affinity weights
 *
 * Parameters:
 * seg_ptr - on input a pointer to the segmentation volume before merging,
 *           on output, the merged segmentation
 * rg_ptr - a pointer to the region graph pairs. This should be the maximal
 *          spanning tree of the original region graph.
 */
template< typename ID, typename F> inline void mergerg(
    const volume_ptr<ID>& seg_ptr,
    const region_graph_ptr<ID, F> rg_ptr,
    F thd) {
  std::map<ID, std::tuple<ID, F>> pd;
  std::cout << "Mergerg()" << std::endl << std::flush;
  size_t num = 0;
  for (auto e:*rg_ptr) {
    ID c = std::get<1>(e);
    ID p = std::get<2>(e);
    F a = std::get<0>(e);
    if (a > thd) {
      ++num;
    }
    pd.emplace(c, std::tuple<ID, F>(p, a));
  }
  std::cout << "Found " << num << " worthy edges." << std::endl << std::flush;
  std::map<ID, ID> rd;
  std::set<ID> rset;
  for (auto e:*rg_ptr) {
    F a = std::get<0>(e);
    ID c0 = std::get<1>(e);
    ID p0 = std::get<2>(e);
    ID p = p0;
    ID c = c0;
    while ((a >= thd) && (pd.count(p) > 0)) {
      a = std::get<1>(pd[p]);
      p = std::get<0>(pd[p]);
    }
    if (p != p0) {
      rd[c0] = p;
      rset.insert(p);
    }
  }
  std::cout << "Finished building rset and rd" << std::endl 
            << "Relabeling voxels" << std::endl << std::flush;
  num = 0;
  for (auto plane:*seg_ptr) {
    for (auto raster:plane) {
      for (auto voxel:raster) {
        if (rd.count(voxel) > 0) {
          voxel = rd[voxel];
          ++num;
        }
      }
    }
  }
  std::cout << "Relabeled " << num << " voxels." << std::endl << std::flush;
}