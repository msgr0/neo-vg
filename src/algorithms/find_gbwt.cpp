/**
 * \file find_gbwt.cpp
 */

#include "find_gbwt.hpp"
#include <vg/io/vpkg.hpp>

#include "../gbzgraph.hpp"

namespace vg {
namespace algorithms {

/// Find the GBWTGraph that is part of the given handle graph (a GBWTGraph or
/// GBZGraph), or null if none exists.
static const gbwtgraph::GBWTGraph* find_gbwtgraph(const HandleGraph* graph) {
    if (!graph) {
        // No graph means no translation.
        return nullptr;
    }
    if (dynamic_cast<const gbwtgraph::GBWTGraph*>(graph)) {
        // If it already is one, return it
        return dynamic_cast<const gbwtgraph::GBWTGraph*>(graph);
    }
    if (dynamic_cast<const GBZGraph*>(graph)) {
        // If it's a GBZGraph, go get the GBWTGraph and return that.
        return &dynamic_cast<const GBZGraph*>(graph)->gbz.graph;
    }
    // Otherwise there's no applicable GBWTGraph
    return nullptr;
}

const gbwt::GBWT* find_gbwt(const HandleGraph* graph) {
    const gbwtgraph::GBWTGraph* typed_graph = find_gbwtgraph(graph);
    if (!typed_graph) {
        return nullptr;
    }
    return typed_graph->index;
}

const gbwt::GBWT* find_gbwt(const HandleGraph* graph, std::unique_ptr<gbwt::GBWT>& holder, const std::string& filename) {
    if (!filename.empty()) {
        holder = vg::io::VPKG::load_one<gbwt::GBWT>(filename);

        if (holder.get() == nullptr) {
          // Complain if we couldn't get it but were supposed to.
          cerr << "error:[vg::algorithms::find_gbwt] unable to load gbwt index file " << filename << endl;
          exit(1);
        }
        
        return holder.get();
    }
    // If we don't need to load it, try and get it from the graph.
    return find_gbwt(graph);
}

}
}
