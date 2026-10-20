#ifndef VG_ALGORITHMS_THREE_EDGE_CONNECTED_COMPONENTS_HPP_INCLUDED
#define VG_ALGORITHMS_THREE_EDGE_CONNECTED_COMPONENTS_HPP_INCLUDED

#include <functional>

#include "../handle.hpp"

namespace vg {
namespace algorithms {

using namespace std;

// Interface

/**
 * Get the three-edge-connected components of an arbitrary graph (not
 * necessarily a handle graph). Only recognizes one kind of edge and one kind
 * of node. Nodes are handles (treated as opaque values).
 *
 * Takes a function that loops an iteratee over all nodes, and a function that,
 * given a node, loops an iteratee over all nodes connected to it.
 *
 * Calls same_component with pairs of nodes in (at least) a spanning tree of
 * the set of nodes in each component (not restricted to the input graph).
 * Doing merge operations on a union-find can get you the set of components.
 * The callback MUST NOT modify the graph!
 *
 * If you have a graph where you can easily rank the nodes, don't use this. Use
 * three_edge_connected_components_dense() instead. The first thing this
 * function does is asign nodes a dense, 0-based rank space.  
 */
void three_edge_connected_component_merges(const function<void(const function<void(handle_t)>&)>& for_each_node,
    const function<void(handle_t, const function<void(handle_t)>&)>& for_each_connected_node,
    const function<void(handle_t, handle_t)>& same_component);


/**
 * Get the three-edge-connected components of an arbitrary graph (not
 * necessarily a handle graph). Only recognizes one kind of edge and one kind
 * of node. Nodes are dense positive integers starting with 0.
 *
 * Takes a total node count, a suggested root (or 0), and a function that,
 * given a node, loops an iteratee over all nodes connected to it.
 *
 * Calls same_component with pairs of nodes in (at least) a spanning tree of
 * the set of nodes in each component (not restricted to the input graph).
 * Doing merge operations on a union-find can get you the set of components.
 * The callback MUST NOT modify the graph!
 */
void three_edge_connected_component_merges_dense(size_t node_count, size_t first_root,
    const function<void(size_t, const function<void(size_t)>&)>& for_each_connected_node,
    const function<void(size_t, size_t)>& same_component);

    
    
}
}

#endif
