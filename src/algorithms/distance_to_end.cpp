#include "distance_to_end.hpp"

#include <unordered_set>

namespace vg {
namespace algorithms {

using namespace std;

/// True if the node has no edges on its left (go_left) or right side.
static bool is_end_node(handle_t h, bool go_left, const HandleGraph* g) {
    bool no_edges = true;
    g->follow_edges(g->forward(h), go_left, [&](const handle_t& ignored) {
        // We only need one
        no_edges = false;
        return false;
    });
    return no_edges;
}

static int32_t distance_to_end(handle_t h, int32_t limit, int32_t dist, bool go_left,
                               unordered_set<handle_t>& seen, const HandleGraph* graph) {
    if (seen.count(h)) return -1;
    seen.insert(h);
    if (limit <= 0) {
        return -1;
    }

    if (is_end_node(h, go_left, graph)) {
        return dist;
    }

    int32_t t = -1;

    graph->follow_edges(h, go_left, [&](const handle_t& current) {
        int32_t l = graph->get_length(current);
        t = distance_to_end(current, limit-l, dist+l, go_left, seen, graph);
        return t == -1;
    });
    return t;
}

int32_t distance_to_end(handle_t h, int32_t limit, bool go_left, const HandleGraph* graph) {
    unordered_set<handle_t> seen;
    return distance_to_end(h, limit, 0, go_left, seen, graph);
}

}
}
