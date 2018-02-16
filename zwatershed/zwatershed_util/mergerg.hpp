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
 * seg - on input a reference to the segmentation volume before merging,
 *           on output, the merged segmentation
 * rg_ptr - a pointer to the region graph pairs. This should be the maximal
 *          spanning tree of the original region graph.
 */
template< typename VOLUME_T, typename F> inline void mergerg(
    VOLUME_T& seg,
    const region_graph_ptr< typename VOLUME_T::element, F> rg_ptr,
    F thd) {
  using ID = typename VOLUME_T::element;
  std::map<ID, std::tuple<ID, F>> pd;
  std::cout << "Mergerg()" << std::endl << std::flush;
  size_t num = 0;
  for (auto e:*rg_ptr) {
    F a = std::get<0>(e);
    ID c = std::get<1>(e);
    ID p = std::get<2>(e);
    if (a > thd) {
      ++num;
    }
    pd.emplace(c, std::tuple<ID, F>(p, a));
  }
  std::cout << "aa-Found " << num << " worthy edges." << std::endl << std::flush;
  std::map<ID, ID> rd;
  std::set<ID> rset;
  for (auto e:*rg_ptr) {
    F a = std::get<0>(e);
    if (a < thd) continue;
    ID c0 = std::get<1>(e);
    ID p0 = std::get<2>(e);
    // std::cout << a << ","<<c0 << "," << p0<<std::endl << std::flush;
    ID p = p0;
    ID c = c0;
    while (pd.count(p) > 0) {
      a = std::get<1>(pd[p]);
      if (a < thd) break;
      p = std::get<0>(pd[p]);
    }
    rd[c0] = p;
    rset.insert(p);
  }
  std::cout << "Found " << rset.size() << " parents" << std::endl 
            << "Merged " << rd.size() << " children" << std::endl
            << "Relabeling voxels" << std::endl << std::flush;
  num = 0;
  for (auto plane:seg) {
    for (auto raster:plane) {
      for (ID &voxel:raster) {
        if (rd.count(voxel) > 0) {
          voxel = rd[voxel];
          ++num;
        }
      }
    }
  }
  std::cout << "Relabeled " << num << " voxels." << std::endl << std::flush;
}
