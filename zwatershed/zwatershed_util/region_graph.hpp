#pragma once

#include "types.hpp"

#include <cstddef>
#include <iostream>
#include <map>
#include <utility>

template< typename VOLUME_T, typename AFF_T >
inline region_graph_ptr<typename VOLUME_T::element, typename AFF_T::element>
get_region_graph( const AFF_T & aff,
                  const VOLUME_T & seg,
                  std::size_t max_segid)
{
    using ID = typename VOLUME_T::element;
    using F = typename AFF_T::element;
    std::ptrdiff_t xdim = aff.shape()[0];
    std::ptrdiff_t ydim = aff.shape()[1];
    std::ptrdiff_t zdim = aff.shape()[2];

    region_graph_ptr<ID,F> rg_ptr( new region_graph<ID,F> );

    region_graph<ID,F>& rg = *rg_ptr;

    std::vector<std::map<ID,F>> edges(max_segid+1);

    for ( std::ptrdiff_t z = 0; z < zdim; ++z )
        for ( std::ptrdiff_t y = 0; y < ydim; ++y )
            for ( std::ptrdiff_t x = 0; x < xdim; ++x )
            {
                if ( (x > 0) && seg[x][y][z] && seg[x-1][y][z] )
                {
                    auto mm = std::minmax(seg[x][y][z], seg[x-1][y][z]);
                    F& curr = edges[mm.first][mm.second];
                    curr = std::max(curr, aff[x][y][z][0]);
                }
                if ( (y > 0) && seg[x][y][z] && seg[x][y-1][z] )
                {
                    auto mm = std::minmax(seg[x][y][z], seg[x][y-1][z]);
                    F& curr = edges[mm.first][mm.second];
                    curr = std::max(curr, aff[x][y][z][1]);
                }
                if ( (z > 0) && seg[x][y][z] && seg[x][y][z-1] )
                {
                    auto mm = std::minmax(seg[x][y][z], seg[x][y][z-1]);
                    F& curr = edges[mm.first][mm.second];
                    curr = std::max(curr, aff[x][y][z][2]);
                }
            }

    for ( ID id1 = 1; id1 <= max_segid; ++id1 )
    {
        for ( const auto& p: edges[id1] )
        {
            rg.emplace_back(p.second, id1, p.first);
            // std::cout << p.second << " " << id1 << " " << p.first << "\n";
        }
    }

    std::cout << "Region graph size: " << rg.size() << std::endl;

    std::stable_sort(std::begin(rg), std::end(rg),
                     std::greater<std::tuple<F,ID,ID>>());

    return rg_ptr;
}

template< typename VOLUME_T, typename AFF_T >
inline region_graph_ptr<typename VOLUME_T::element, typename AFF_T::element>
get_region_graph_average( const AFF_T & aff,
                          const VOLUME_T & seg,
                          std::size_t max_segid)
{
    using ID = typename VOLUME_T::element;
    using F = typename AFF_T::element;
    std::ptrdiff_t xdim = aff.shape()[0];
    std::ptrdiff_t ydim = aff.shape()[1];
    std::ptrdiff_t zdim = aff.shape()[2];
    region_graph_ptr<ID,F> rg_ptr( new region_graph<ID,F> );

    region_graph<ID,F>& rg = *rg_ptr;

    using EDGE_T = std::pair<F, std::uint64_t>;
    std::vector<std::map<ID,EDGE_T> > edges(max_segid+1);

    for ( std::ptrdiff_t z = 0; z < zdim; ++z )
        for ( std::ptrdiff_t y = 0; y < ydim; ++y )
            for ( std::ptrdiff_t x = 0; x < xdim; ++x )
            {
                if ( (x > 0) && seg[x][y][z] && seg[x-1][y][z] )
                {
                    auto mm = std::minmax(seg[x][y][z], seg[x-1][y][z]);
                    EDGE_T& curr = edges[mm.first][mm.second];
                    curr.first += aff[x][y][z][0];
                    ++curr.second;
                }
                if ( (y > 0) && seg[x][y][z] && seg[x][y-1][z] )
                {
                    auto mm = std::minmax(seg[x][y][z], seg[x][y-1][z]);
                    EDGE_T& curr = edges[mm.first][mm.second];
                    curr.first += aff[x][y][z][1];
                    ++curr.second;
                }
                if ( (z > 0) && seg[x][y][z] && seg[x][y][z-1] )
                {
                    auto mm = std::minmax(seg[x][y][z], seg[x][y][z-1]);
                    EDGE_T& curr = edges[mm.first][mm.second];
                    curr.first += aff[x][y][z][2];
                    ++curr.second;
                }
            }

    for ( ID id1 = 1; id1 <= max_segid; ++id1 )
    {
        for ( const auto& p: edges[id1] )
        {
            F avg = p.second.first / ((float) p.second.second);
            rg.emplace_back(avg, id1, p.first);
            // std::cout << p.second << " " << id1 << " " << p.first << "\n";
        }
    }

    std::cout << "Region graph size: " << rg.size() << std::endl;

    std::stable_sort(std::begin(rg), std::end(rg),
                     std::greater<std::tuple<F,ID,ID>>());

    return rg_ptr;
}
