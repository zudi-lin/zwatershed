#pragma once

#include "types.hpp"
#include "mst.hpp"
#include "mergerg.hpp"
#include <zi/disjoint_sets/disjoint_sets.hpp>
#include <map>
#include <vector>
#include <set>
template< typename VOLUME_T, typename F, typename M >
inline region_graph_ptr<typename VOLUME_T::element, F> 
merge_segments_with_function_dw( 
    VOLUME_T &seg,
    const region_graph<typename VOLUME_T::element,F> &rg,
    std::vector<std::size_t>& counts,
    const M& size_th,
    const F& weight_th,
    const M& lowt,
    const F& merge_th,
    const M do_mst)
{
    using ID=typename VOLUME_T::element;
    zi::disjoint_sets<ID> sets(counts.size());

    for ( auto& it: rg )
    {
        F weight = std::get<0>(it);

        ID s1 = sets.find_set(std::get<1>(it));
        ID s2 = sets.find_set(std::get<2>(it));

        // std::cout << s1 << " " << s2 << " " << size << "\n";

        if ( s1 != s2 && s1 && s2 && weight>weight_th )
        {
            if ( (counts[s1] < size_th) || (counts[s2] < size_th) )
            {
                counts[s1] += counts[s2];
                counts[s2]  = 0;
                ID s = sets.join(s1,s2);
                std::swap(counts[s], counts[s1]);
            }
        }
    }
    
    // counts: start from 0 
    std::vector<ID> remaps(counts.size(), 0);

    ID next_id = 1;

    std::size_t low = static_cast<std::size_t>(lowt);

    for ( ID id = 0; id < counts.size(); ++id )
    {
        ID s = sets.find_set(id);
        if ( s && (remaps[s] == 0) && (counts[s] >= low) )
        {
            remaps[s] = next_id;
            counts[next_id] = counts[s];
            ++next_id;
        }
    }
    std::cout<<"next_id: "<<counts.size()<<","<<next_id<<std::endl;

    counts.resize(next_id);

    std::ptrdiff_t xdim = seg.shape()[0];
    std::ptrdiff_t ydim = seg.shape()[1];
    std::ptrdiff_t zdim = seg.shape()[2];

    std::ptrdiff_t total = xdim * ydim * zdim;

    for (auto plane:seg) {
      for (auto raster:plane) {
        for (ID &voxel:raster) {
            voxel = remaps[sets.find_set(voxel)];
        }
      }
    }

    std::cout << "\tDone with remapping, total: " << (next_id-1) << std::endl;

    region_graph_ptr<ID,F> new_rg_ptr(new region_graph<ID, F>);
    region_graph<ID, F> &new_rg=*new_rg_ptr;

    std::vector<std::set<ID>> in_rg(next_id);

    for ( auto& it: rg )
    {
        ID s1 = remaps[sets.find_set(std::get<1>(it))];
        ID s2 = remaps[sets.find_set(std::get<2>(it))];

        if ( s1 != s2 && s1 && s2 )
        {
            auto mm = std::minmax(s1,s2);
            if ( in_rg[mm.first].count(mm.second) == 0 )
            {
                // std::cout<<sets.find_set(std::get<1>(it))<<","<<sets.find_set(std::get<2>(it))<<","<<s1<<","<<s2<<std::endl;
                new_rg.push_back(std::make_tuple(std::get<0>(it), mm.first, mm.second));
                in_rg[mm.first].insert(mm.second);
            }
        }
    }
    if (do_mst>0){
        std::cout << "Do MST" << std::endl;
        new_rg_ptr = mst(new_rg_ptr, counts.size());
        std::cout << "New region graph size after mst = " << new_rg_ptr->size() << std::endl;
        mergerg(seg, new_rg_ptr, merge_th);
        std::cout << "New region graph size after mergerg = " << new_rg_ptr->size() << std::endl;

        std::cout << "\tDone with updating the region graph, size: "
                  << rg.size() << std::endl;
    }
    return new_rg_ptr;
}

