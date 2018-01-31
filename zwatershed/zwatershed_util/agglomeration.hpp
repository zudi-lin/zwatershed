#pragma once

#include "types.hpp"
#include "mst.hpp"
#include "mergerg.hpp"
#include <zi/disjoint_sets/disjoint_sets.hpp>
#include <map>
#include <vector>
#include <set>
template< typename ID, typename F, typename M >
inline region_graph_ptr<ID, F> merge_segments_with_function_dw( const volume_ptr<ID>& seg_ptr,
                                          const region_graph_ptr<ID,F> rg_ptr,
                                          std::vector<std::size_t>& counts,
                                          const M& size_th,
                                          const F& weight_th,
                                          const M& lowt,
                                          const F& merge_th)
{
    zi::disjoint_sets<ID> sets(counts.size());

    region_graph<ID,F>& rg  = *rg_ptr;

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

    std::vector<ID> remaps(counts.size());

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

    counts.resize(next_id);

    std::ptrdiff_t xdim = seg_ptr->shape()[0];
    std::ptrdiff_t ydim = seg_ptr->shape()[1];
    std::ptrdiff_t zdim = seg_ptr->shape()[2];

    std::ptrdiff_t total = xdim * ydim * zdim;

    ID* seg_raw = seg_ptr->data();

    for ( std::ptrdiff_t idx = 0; idx < total; ++idx )
    {
        seg_raw[idx] = remaps[sets.find_set(seg_raw[idx])];
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
                new_rg.push_back(std::make_tuple(std::get<0>(it), mm.first, mm.second));
                in_rg[mm.first].insert(mm.second);
            }
        }
    }

    new_rg_ptr = mst(new_rg_ptr, counts.size());
    std::cout << "New region graph size after mst = " << new_rg_ptr->size() << std::endl;
    mergerg(seg_ptr, new_rg_ptr, merge_th);
    std::cout << "New region graph size after mergerg = " << new_rg_ptr->size() << std::endl;

    std::cout << "\tDone with updating the region graph, size: "
              << rg.size() << std::endl;
    return new_rg_ptr;
}
template< typename ID, typename F, typename FN, typename M >
inline void merge_segments_with_function( const volume_ptr<ID>& seg_ptr,
                                          const region_graph_ptr<ID,F> rg_ptr,
                                          std::vector<std::size_t>& counts,
                                          const FN& func,
                                          const M& lowt,
                                           bool recreate_rg)
{
    zi::disjoint_sets<ID> sets(counts.size());

    region_graph<ID,F>& rg  = *rg_ptr;

    for ( auto& it: rg )
    {
        std::size_t size = func(std::get<0>(it));

        if ( size == 0 )
        {
            break;
        }

        ID s1 = sets.find_set(std::get<1>(it));
        ID s2 = sets.find_set(std::get<2>(it));

        // std::cout << s1 << " " << s2 << " " << size << "\n";

        if ( s1 != s2 && s1 && s2 )
        {
            if ( (counts[s1] < size) || (counts[s2] < size) )
            {
                counts[s1] += counts[s2];
                counts[s2]  = 0;
                ID s = sets.join(s1,s2);
                std::swap(counts[s], counts[s1]);
            }
        }
    }

    std::vector<ID> remaps(counts.size());

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

    counts.resize(next_id);

    std::ptrdiff_t xdim = seg_ptr->shape()[0];
    std::ptrdiff_t ydim = seg_ptr->shape()[1];
    std::ptrdiff_t zdim = seg_ptr->shape()[2];

    std::ptrdiff_t total = xdim * ydim * zdim;

    ID* seg_raw = seg_ptr->data();

    for ( std::ptrdiff_t idx = 0; idx < total; ++idx )
    {
        seg_raw[idx] = remaps[sets.find_set(seg_raw[idx])];
    }

    std::cout << "\tDone with remapping, total: " << (next_id-1) << std::endl;

    region_graph<ID,F> new_rg;

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
                new_rg.push_back(std::make_tuple(std::get<0>(it), mm.first, mm.second));
                in_rg[mm.first].insert(mm.second);
            }
        }
    }
    if(recreate_rg)
        rg.swap(new_rg);

    std::cout << "\tDone with updating the region graph, size: "
              << rg.size() << std::endl;
}
