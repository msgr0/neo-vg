#ifndef VG_HAPLOTYPE_EXTRACTER_HPP_INCLUDED
#define VG_HAPLOTYPE_EXTRACTER_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <gbwt/gbwt.h>
#include <vg/vg.pb.h>
#include "bdsg/hash_graph.hpp"

#include "algorithms/expand_context.hpp"
#include "gbwt_helper.hpp"

namespace vg {

using namespace std;

using thread_t = vector<gbwt::node_type>;
    
// Walk forward from a node, collecting all haplotypes.  Also do a regular
// subgraph search for all the paths too.  Haplotype thread i will be embedded
// as Paths a path with name thread_i.  Each path name (including threads) is
// mapped to a frequency in out_thread_frequencies.  Haplotypes will be pulled
// from the GBWT index.
void trace_haplotypes_and_paths(const PathHandleGraph& source,
                                const gbwt::GBWT& haplotype_database,
                                vg::id_t start_node, int extend_distance,
                                Graph& out_graph,
                                map<string, int>& out_thread_frequencies,
                                bool expand_graph = true);

// Turns a (GBWT-based) thread_t into a (vg-based) Path
Path path_from_thread_t(thread_t& t, const HandleGraph& source);

// Lists all the sub-haplotypes of nodes starting at
// node start_node from the set of haplotypes embedded in the geven GBWT
// haplotype database.  At each step stop_fn() is called on the thread being created, and if it returns true
// then the search stops and the thread is added two the list to be returned.
vector<pair<vector<gbwt::node_type>, gbwt::SearchState> > list_haplotypes(const HandleGraph& graph,
                                                                          const gbwt::GBWT& gbwt,
                                                                          handle_t start,
                                                                          function<bool(const vector<gbwt::node_type>&)> stop_fn);

}

#endif
