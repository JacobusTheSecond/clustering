/**
 * @file shortest_path.cc
 *
 * Implement the BGL shortest path algorithm wrappers.
 */

/*
 * David Gleich
 * 19 April 2006
 */

/*
 * 18 April 2007
 * Added src/dst vertex pairs to all the calls to allow partial searches.
 * Corrected small documentation bugs.
 *
 * 9 July 2007
 * Switched to simple_csr_matrix graph type
 */

#include "include/matlab_bgl.h"

#include "yasmic/simple_csr_matrix_as_graph.hpp"
#include "yasmic/iterator_utility.hpp"

#include <boost/graph/dag_shortest_paths.hpp>

#include "visitor_macros.hpp"
#include "stop_visitors.hpp"
#include "libmbgl_util.hpp"

struct stop_dag {}; // stop dag exception

int dag_sp(
    mbglIndex nverts, mbglIndex *ja, mbglIndex *ia, double *weight, /* connectivity params */
    mbglIndex src, mbglIndex dst, /* problem data */
    double* d, mbglIndex *pred, double dinf)
{
    using namespace yasmic;
    using namespace boost;

    typedef simple_csr_matrix<mbglIndex,double> crs_weighted_graph;
    crs_weighted_graph g(nverts, nverts, ia[nverts], ia, ja, weight);

    if (dst == nverts) {
        dag_shortest_paths(g, src,
            distance_inf(dinf).predecessor_map(pred).distance_map(d));
    } else {
        try {
            dag_shortest_paths(g, src,
                distance_inf(dinf).predecessor_map(pred).distance_map(d).
                visitor(make_dijkstra_visitor(
                    stop_search_on_vertex_target(dst, stop_dag(), on_discover_vertex()))));
        } catch (stop_dag) {}
    }

    return (0);
}
