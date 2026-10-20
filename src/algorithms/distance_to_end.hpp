#ifndef VG_ALGORITHMS_DISTANCE_TO_END_HPP_INCLUDED
#define VG_ALGORITHMS_DISTANCE_TO_END_HPP_INCLUDED

#include "../handle.hpp"

namespace vg {
namespace algorithms {

using namespace std;

/// Get the distance in bases from the node to the closest head node (go_left =
/// true) or tail node (go_left = false) of the graph, counting the bases of
/// each node stepped onto, or -1 if that distance reaches the limit or there is
/// no such node.
int32_t distance_to_end(handle_t h, int32_t limit, bool go_left, const HandleGraph* graph);

}
}

#endif
