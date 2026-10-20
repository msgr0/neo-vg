#include <iostream>
#include "haplotype_extracter.hpp"
#include "graph.hpp"

namespace vg {

using namespace std;

void trace_haplotypes_and_paths(const PathHandleGraph& source, const gbwt::GBWT& haplotype_database,
                                vg::id_t start_node, int extend_distance,
                                Graph& out_graph,
                                map<string, int>& out_thread_frequencies,
                                bool expand_graph) {
  // get our haplotypes
  handle_t n = source.get_handle(start_node, false);
  vector<pair<thread_t, gbwt::SearchState> > haplotypes = list_haplotypes(source, haplotype_database, n,
                                                                          [&extend_distance](const vector<gbwt::node_type>& new_thread) {
                                                                              return new_thread.size() >= extend_distance;
                                                                          });

#ifdef debug
  cerr << "Haplotype database " << &haplotype_database << " produced " << haplotypes.size() << " haplotypes" << endl;
#endif

  if (expand_graph) {
      // get our subgraph and "regular" paths by expanding forward
      handle_t handle = source.get_handle(start_node);
      bdsg::HashGraph extractor;
      extractor.create_handle(source.get_sequence(handle), source.get_id(handle));
      // TODO: is expanding only forward really the right behavior here?
      algorithms::expand_context_with_paths(&source, &extractor, extend_distance, true, true, false);
      
      // Convert to Protobuf I guess
      from_path_handle_graph(extractor, out_graph);
  }

  // add a frequency of 1 for each normal path
  for (int i = 0; i < out_graph.path_size(); ++i) {
    out_thread_frequencies[out_graph.path(i).name()] = 1;
  }

  // add our haplotypes to the subgraph, naming ith haplotype "thread_i"
  for (int i = 0; i < haplotypes.size(); ++i) {
    Path p = path_from_thread_t(haplotypes[i].first, source);
    p.set_name("thread_" + to_string(i));
    out_thread_frequencies[p.name()] = haplotypes[i].second.size();
    *(out_graph.add_path()) = std::move(p);
  }
}

Path path_from_thread_t(thread_t& t, const HandleGraph& source) {
	Path toReturn;
	int rank = 1;
	for(int i = 0; i < t.size(); i++) {
		Mapping* mapping = toReturn.add_mapping();

        // Set up the position
        mapping->mutable_position()->set_node_id(gbwt::Node::id(t[i]));
        mapping->mutable_position()->set_is_reverse(gbwt::Node::is_reverse(t[i]));

        // Set up the edits
        Edit* e = mapping->add_edit();
        size_t l = source.get_length(source.get_handle(gbwt::Node::id(t[i])));
        e->set_from_length(l);
        e->set_to_length(l);

        // Set the rank
        mapping->set_rank(rank++);
    }
    // We're done making the path
    return toReturn;
}

vector<pair<vector<gbwt::node_type>, gbwt::SearchState> > list_haplotypes(const HandleGraph& graph,
                                                                          const gbwt::GBWT& gbwt,
                                                                          handle_t start,
                                                                          function<bool(const vector<gbwt::node_type>&)> stop_fn) {
    
    // Keep track of all the different paths we're extending
    vector<pair<vector<gbwt::node_type>, gbwt::SearchState> > search_intermediates;
    vector<pair<vector<gbwt::node_type>, gbwt::SearchState> > search_results;

    // Look up the start node in GBWT and start a thread
    gbwt::node_type start_node = handle_to_gbwt(graph, start);    
    vector<gbwt::node_type> first_thread = {start_node};
    gbwt::SearchState first_state = gbwt.find(start_node);
    
#ifdef debug
    cerr << "Start with state " << first_state << " for node " << gbwt::Node::id(start_node)  << ":"
         << gbwt::Node::is_reverse(start_node) << endl;
#endif

    if (!first_state.empty()) {
        search_intermediates.push_back(make_pair(first_thread, first_state));
    }

    while(!search_intermediates.empty()) {

        // pick up a thread to continue from the queue
        auto last = std::move(search_intermediates.back());
        search_intermediates.pop_back();

        vector<tuple<handle_t, gbwt::node_type, gbwt::SearchState>> next_handle_states;
        graph.follow_edges(gbwt_to_handle(graph, last.first.back()), false, [&](const handle_t& next) {
                // extend the last node of the thread using gbwt
                auto extend_node = handle_to_gbwt(graph, next);
                auto new_state = gbwt.extend(last.second, extend_node);
#ifdef debug
                cerr << "Extend state " << last.second << " to " << new_state << " with " << gbwt::Node::id(extend_node) << endl;
#endif
                if (!new_state.empty()) {
                    next_handle_states.push_back(make_tuple(next, extend_node, new_state));
                }                    
            });

        for (auto& nhs : next_handle_states) {
            
            const handle_t& next = get<0>(nhs);
            gbwt::node_type& extend_node = get<1>(nhs);
            gbwt::SearchState& new_state = get<2>(nhs);
                
            vector<gbwt::node_type> new_thread;
            if (&nhs == &next_handle_states.back()) {
                // avoid a copy by re-using the vector for the last thread. this way simple cases
                // like scanning along one path don't blow up to n^2
                new_thread = std::move(last.first);
            } else {
                new_thread = last.first;
            }                        
            new_thread.push_back(extend_node);

            if (stop_fn(new_thread)) {
#ifdef debug
                cerr << "\tGot " << new_state.size() << " results at limit; emitting" << endl;
#endif
                search_results.push_back(make_pair(std::move(new_thread), new_state));
            }
            else {
#ifdef debug
                cerr << "\tGot " << new_state.size() << " results; extending more" << endl;
#endif
                search_intermediates.push_back(make_pair(std::move(new_thread), new_state));
            }
        }
    }
    
    return search_results;
}



}
