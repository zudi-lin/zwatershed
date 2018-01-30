#pragma once

#include "types.hpp"

#include <cstddef>
#include <iostream>
#include <map>
using namespace std;
template< typename ID, typename F >
inline region_graph_ptr<ID,F>
get_region_graph_arb( //const affinity_graph_ptr<F>& aff_ptr,
                const ID*node1, const ID*node2, const F*edgeWeight, int n_edge,
                  const volume_ptr<ID> seg_ptr, std::size_t max_segid)
{

    region_graph_ptr<ID,F> rg_ptr( new region_graph<ID,F> );
    region_graph<ID,F>& rg = *rg_ptr;
    ID* seg = seg_ptr->data();
    std::vector<std::map<ID,F>> edges(max_segid+1);
    std::ptrdiff_t x = 1;
    for(int i=0;i<n_edge;i++){
        ID n1 = node1[i];
        ID n2 = node2[i];
        F w = edgeWeight[i];
        auto mm = std::minmax(seg[n1],seg[n2]); //need to make sure this is in right order
        F& curr = edges[mm.first][mm.second];
        curr = std::max(curr, w);
    }


    for ( ID id1 = 1; id1 <= max_segid; ++id1 )
        for ( const auto& p: edges[id1] )
            rg.emplace_back(p.second, id1, p.first);


    std::stable_sort(std::begin(rg), std::end(rg), std::greater<std::tuple<F,ID,ID>>());


    return rg_ptr;
}

