#include "traversal_finder.hpp"
#include "genotypekit.hpp"
#include "algorithms/k_widest_paths.hpp"
#include "cactus.hpp"
#include "gbwt_helper.hpp"
#include "haplotype_extracter.hpp"
//#define debug

namespace vg {

using namespace std;

string traversal_to_string(const PathHandleGraph* graph, const Traversal& traversal, int64_t max_steps) {
    string s;
    function<string(handle_t)> handle_to_string = [&](handle_t handle) {
        string ss = graph->get_is_reverse(handle) ? "<" : ">";
        return ss + std::to_string(graph->get_id(handle));
    };
    if (traversal.size() <= max_steps) {
        for (const handle_t& handle : traversal) {
            s += handle_to_string(handle);
        }
    } else {
        for (int64_t i = 0; i < max_steps / 2; ++i) {
            s += handle_to_string(traversal[i]);
        }
        s += "...";
        for (int64_t i = 0; i < max_steps / 2; ++i) {
            s += handle_to_string(traversal[traversal.size() - 1 - i]);
        }
    }
    return s;
}

string graph_interval_to_string(const HandleGraph* graph, const handle_t& start_handle, const handle_t& end_handle) {
    function<string(handle_t)> handle_to_string = [&](handle_t handle) {
        string ss = graph->get_is_reverse(handle) ? "<" : ">";
        return ss + std::to_string(graph->get_id(handle));
    };
    return handle_to_string(start_handle) + handle_to_string(end_handle);
}

string path_interval_to_string(const PathHandleGraph* graph, const PathInterval& path_interval) {
    string s = graph->get_path_name(graph->get_path_handle_of_step(path_interval.first));
    return s + ":" + graph_interval_to_string(graph, graph->get_handle_of_step(path_interval.first),
                                              graph->get_handle_of_step(path_interval.second));   
}

ExhaustiveTraversalFinder::ExhaustiveTraversalFinder(const HandleGraph& graph, SnarlManager& snarl_manager,
                                                     bool include_reversing_traversals) :
    graph(graph), snarl_manager(snarl_manager),
    include_reversing_traversals(include_reversing_traversals) {
    // nothing more to do
}
    
ExhaustiveTraversalFinder::~ExhaustiveTraversalFinder() {
    // no heap objects
}

void ExhaustiveTraversalFinder::stack_up_valid_walks(handle_t walk_head, vector<Visit>& stack) {
    
    id_t head_id = graph.get_id(walk_head);

    // get all edges involving this node so we can filter them down to valid walks
    graph.follow_edges(walk_head, false, [&](const handle_t next_node) {
            if (visit_next_node(next_node)) {
                Visit next_visit;
                next_visit.set_node_id(graph.get_id(next_node));
                next_visit.set_backward(graph.get_is_reverse(next_node));
                stack.push_back(next_visit);
            }
        });
}

void ExhaustiveTraversalFinder::add_traversals(vector<SnarlTraversal>& traversals,
                                               handle_t traversal_start,
                                               unordered_set<handle_t>& stop_at,
                                               unordered_set<handle_t>& yield_at) {
    // keeps track of the walk of the DFS traversal
    list<Visit> path;
    
    // these mark the start of the edges out of the node that is on the head of the path
    // they can be used to see how many nodes we need to peel off the path when we're
    // backtracking
    Visit stack_sentinel;
    
    // initialize stack for DFS traversal of site
    vector<Visit> stack{to_visit(graph, traversal_start)};
    
    while (stack.size()) {
        
        Visit node_visit = stack.back();
        stack.pop_back();
        
        // we have traversed all of edges out of the head of the path, so we can pop it off
        if (node_visit == stack_sentinel) {
            path.pop_back();
            continue;
        }
        
        // have we finished a traversal through the site?
        handle_t node_handle = graph.get_handle(node_visit.node_id(), node_visit.backward());
        if (stop_at.count(node_handle)) {
            if (yield_at.count(node_handle)) {
                // yield path as a snarl traversal
                traversals.emplace_back();
                
                // record the traversal in the return value
                for (auto iter = path.begin(); iter != path.end(); iter++) {
                    *traversals.back().add_visit() = *iter;
                }
                // add the final visit
                *traversals.back().add_visit() = node_visit;
            }
            
            // don't proceed to add more onto the DFS stack
            continue;
        }
        
        // mark the beginning of this node's edges forward in the stack
        stack.push_back(stack_sentinel);
        
        // make a visit through the node traversal and add it to the path
        path.emplace_back();
        path.back().set_node_id(node_visit.node_id());
        path.back().set_backward(node_visit.backward());
        
        // does this traversal point into a child snarl?
        const Snarl* into_snarl = snarl_manager.into_which_snarl(node_visit.node_id(),
                                                                 node_visit.backward());
                                                                 
#ifdef debug
        cerr << "Traversal " << node_visit.node_id() << " " << node_visit.backward() << " enters";
        if (into_snarl != nullptr) {
            cerr << " " << pb2json(*into_snarl) << endl;
        } else {
            cerr << " NULL" << endl; 
        }
#endif
                                                                 
        if (into_snarl && !(node_handle == traversal_start)) {
            // add a visit for the child snarl
            path.emplace_back();
            *path.back().mutable_snarl()->mutable_start() = into_snarl->start();
            *path.back().mutable_snarl()->mutable_end() = into_snarl->end();
            
            // mark the beginning of this child snarls edges forward in the stack
            stack.push_back(stack_sentinel);
            
            // which side of the snarl does the traversal point into?
            if (into_snarl->start().node_id() == node_visit.node_id()
                && into_snarl->start().backward() == node_visit.backward()) {
                // Into the start
#ifdef debug
                cerr << "Entered child through its start" << endl;
#endif
                if (into_snarl->start_end_reachable()) {
                    // skip to the other side and proceed in the orientation that the end node takes.
                    stack.push_back(into_snarl->end());
#ifdef debug
                    cerr << "Stack up " << stack.back().node_id() << " " << stack.back().backward() << endl;
#endif
                }
                
                // if the same side is also reachable, add it to the stack too
                if (into_snarl->start_self_reachable()) {
                    // Make sure to flip it around so we come out of the snarl instead of going in again,
                    Visit rev_visit = into_snarl->start();
                    rev_visit.set_backward(!rev_visit.backward());
                    stack.push_back(rev_visit);
#ifdef debug
                    cerr << "Stack up " << stack.back().node_id() << " " << stack.back().backward() << endl;
#endif
                }
                
            }
            else {
                // Into the end
#ifdef debug
                cerr << "Entered child through its end" << endl;
#endif
                if (into_snarl->start_end_reachable()) {
                    // skip to the other side and proceed in the orientation
                    // *opposite* what the start node takes (i.e. out of the
                    // snarl)
                    Visit rev_visit = into_snarl->start();
                    rev_visit.set_backward(!rev_visit.backward());
                    stack.push_back(rev_visit);
#ifdef debug
                    cerr << "Stack up " << stack.back().node_id() << " " << stack.back().backward() << endl;
#endif
                }
                
                // if the same side is also reachable, add it to the stack too
                if (into_snarl->end_self_reachable()) {
                    Visit rev_visit = into_snarl->end();
                    rev_visit.set_backward(!rev_visit.backward());
                    stack.push_back(rev_visit);
#ifdef debug
                    cerr << "Stack up " << stack.back().node_id() << " " << stack.back().backward() << endl;
#endif
                }
            }
        }
        else {
            // add all of the node traversals we can reach through valid walks to stack
            stack_up_valid_walks(node_handle, stack);
        }
    }
}
    
vector<SnarlTraversal> ExhaustiveTraversalFinder::find_traversals(const Snarl& site) {

    vector<SnarlTraversal> to_return;

    handle_t site_end = graph.get_handle(site.end().node_id(), site.end().backward());
    handle_t site_start = graph.get_handle(site.start().node_id(), site.start().backward());
    handle_t site_rev_start = graph.get_handle(site.start().node_id(), !site.start().backward());
    
    // stop searching when the traversal is leaving the site
    unordered_set<handle_t> stop_at;
    stop_at.insert(site_end);
    stop_at.insert(site_rev_start);
    
    // choose which side(s) can be the end of the traversal
    unordered_set<handle_t> yield_at;
    yield_at.insert(site_end);
    if (include_reversing_traversals) {
        yield_at.insert(site_rev_start);
    }
    
    // search forward from the start and add any traversals that leave the indicated boundaries
    add_traversals(to_return, site_start, stop_at, yield_at);

    if (site.end_self_reachable() && include_reversing_traversals) {
        // if the end is reachable from itself, also look for traversals that both enter and
        // leave through the end
        yield_at.erase(site_rev_start);
        add_traversals(to_return, graph.get_handle(graph.get_id(site_end), !graph.get_is_reverse(site_end)),
                       stop_at, yield_at);
    }
    
    return to_return;
}

PathTraversalFinder::PathTraversalFinder(const PathHandleGraph& graph,
                                         const vector<string>& path_names)  :
    graph(graph) {
    for (const string& path_name : path_names) {
        assert(graph.has_path(path_name));
        paths.insert(graph.get_path_handle(path_name));
    }
}

vector<SnarlTraversal> PathTraversalFinder::find_traversals(const Snarl& site) {
    return find_path_traversals(site).first;
}

vector<Traversal> PathTraversalFinder::find_traversals(const handle_t& snarl_start, const handle_t& snarl_end) {
    return find_path_traversals(snarl_start, snarl_end).first;
}

pair<vector<SnarlTraversal>, vector<pair<step_handle_t, step_handle_t> > > find_path_traversals(const Snarl& site);

pair<vector<SnarlTraversal>, vector<PathInterval>> PathTraversalFinder::find_path_traversals(const Snarl& site) {

    pair<vector<Traversal>, vector<PathInterval>> path_travs = this->find_path_traversals(
        graph.get_handle(site.start().node_id(), site.start().backward()),
        graph.get_handle(site.end().node_id(), site.end().backward()));

    vector<SnarlTraversal> travs(path_travs.first.size());

    for (int64_t i = 0; i < path_travs.first.size(); ++i) {
        SnarlTraversal& proto_trav = travs[i];
        const Traversal& handle_trav = path_travs.first[i];
        for (int64_t j = 0; j < handle_trav.size(); ++j) {
            Visit* visit = proto_trav.add_visit();
            visit->set_node_id(graph.get_id(handle_trav[j]));
            visit->set_backward(graph.get_is_reverse(handle_trav[j]));
        }
    }
    return make_pair(travs, path_travs.second);    
}

pair<vector<Traversal>, vector<PathInterval>> PathTraversalFinder::find_path_traversals(const handle_t& snarl_start, const handle_t& snarl_end) {

    vector<step_handle_t> start_steps = graph.steps_of_handle(snarl_start);
    // try to make this stable across different systems
    std::sort(start_steps.begin(), start_steps.end());
    vector<step_handle_t> end_steps = graph.steps_of_handle(snarl_end);
    
    // use this to skip paths that don't reach the end node
    unordered_set<path_handle_t> end_path_handles;
    for (const step_handle_t& step : end_steps) {
        end_path_handles.insert(graph.get_path_handle_of_step(step));
    }

#ifdef debug
    cerr << "Finding traversals of " << (graph.get_is_reverse(snarl_start) ? "<" : ">") << graph.get_id(snarl_start)
         << (graph.get_is_reverse(snarl_end) ? "<" : ">") << graph.get_id(snarl_end) << " using PathTraversalFinder" << endl
         << " - there are " << start_steps.size() << " start_steps, " << end_steps.size() << " end_steps"
         << " and " << end_path_handles.size() << " end_path_handles" << endl;
#endif

    vector<Traversal> out_travs;
    vector<pair<step_handle_t, step_handle_t> > out_steps;

    // collect nodes on edges that lead out the snarl to make sure we don't search in the wrong direction
    unordered_set<nid_t> block_nodes;
    graph.follow_edges(snarl_start, true, [&](handle_t prev) {
        block_nodes.insert(graph.get_id(prev));
    });
    graph.follow_edges(snarl_end, false, [&](handle_t next) {
        block_nodes.insert(graph.get_id(next));
    });

    for (const step_handle_t& start_step : start_steps) {
        path_handle_t start_path_handle = graph.get_path_handle_of_step(start_step);
        // only crawl paths that have a chance of reaching the end
        if ((paths.empty() || paths.count(start_path_handle)) && end_path_handles.count(start_path_handle)) {

            handle_t end_check = snarl_end;

#ifdef debug
            cerr << " - considering path " << graph.get_path_name(start_path_handle) << endl;
#endif
            // try to make a traversal by walking forward
            Traversal trav;
            bool can_continue = true;
            step_handle_t step = start_step;
            while (can_continue) {
                handle_t handle = graph.get_handle_of_step(step);
                trav.push_back(handle);

                can_continue = false;
                if (graph.has_next_step(step) && handle != snarl_end) {
                    step_handle_t next_step = graph.get_next_step(step);
                    handle_t next_handle = graph.get_handle_of_step(next_step);
                    // todo: we only need this check once
                    if (next_handle == end_check ||
                        !block_nodes.count(graph.get_id(graph.get_handle_of_step(next_step)))) {
                        step = next_step;
                        can_continue = true;
                    }
                }
            }

            if (graph.get_handle_of_step(step) != end_check) {
#ifdef debug
                cerr << "     - failed to find forward traversal of path " << graph.get_path_name(start_path_handle) << endl;
#endif
                // try to make a traversal by walking backward
                end_check = graph.flip(snarl_end);
                
                trav.clear();
                can_continue = true;
                step = start_step;
                while (can_continue) {
                    handle_t handle = graph.flip(graph.get_handle_of_step(step));
                    trav.push_back(handle);

                    can_continue = false;
                    if (graph.has_previous_step(step) && handle != snarl_end) {
                        step_handle_t prev_step = graph.get_previous_step(step);
                        handle_t prev_handle = graph.get_handle_of_step(prev_step);
                        // todo: we only need this check once
                        if (prev_handle == end_check ||
                            !block_nodes.count(graph.get_id(prev_handle))) {
                            step = prev_step;
                            can_continue = true;
                        } 
                    }
                }
            }
            if (graph.get_handle_of_step(step) == end_check) {
                out_travs.push_back(trav);
                out_steps.push_back(make_pair(start_step, step));
            } else {
#ifdef debug
                cerr << "     - failed to find reverse traversal of path " << graph.get_path_name(start_path_handle) << endl;
#endif
                
            }
        }
    }
    
    return make_pair(out_travs, out_steps);
}

TrivialTraversalFinder::TrivialTraversalFinder(const HandleGraph& graph) : graph(graph) {
    // Nothing to do!
}

vector<SnarlTraversal> TrivialTraversalFinder::find_traversals(const Snarl& site) {
    assert(site.start_end_reachable());
    assert(site.directed_acyclic_net_graph());
    
    // We'll fill this in and send it back
    vector<SnarlTraversal> to_return;
    
    // We don't want to be duplicating partial paths, so we store for each
    // NodeTraversal we can reach the previous NodeTraversal we can reach it
    // from.
    unordered_map<handle_t, handle_t> previous;
    
    list<handle_t> stack{graph.get_handle(site.start().node_id(), site.start().backward())};
    
    while (!stack.empty()) { 
        // While there's still stuff on the stack
        
        // Grab the first thing
        handle_t here = stack.front();
        stack.pop_front();
        
        if (graph.get_id(here) == site.end().node_id()) {
            // Trace back a path
            list<handle_t> path;
            
            while (true) {
                // Until we get to the start of the site
                
                // Put this traversal on the front of the path
                path.push_front(here);
                
                if (graph.get_id(here) == site.start().node_id()) {
                    // Stop when we've reached the start of the site
                    break;
                }
                
                // Trace back
                here = previous.at(here);
            }
            
            // Initialize a SnarlTraversal in the return value
            to_return.emplace_back();
            
            // Translate the path into the traversal
            for (handle_t node_traversal : path) {
                *(to_return.back().add_visit()) = to_visit(graph, node_traversal);
            }
            
            // Stop early after having found one path
            break;
        } else {
            // We haven't reached the end of the site

            graph.follow_edges(here, false, [&] (const handle_t& next) {
                    // Look at all the places we can go from this node
                    if (!previous.count(next)) {
                        // Remember how we got there
                        previous[next] = here;
                        // Explore it, depth first
                        stack.push_front(next);
                    }
                });
        }
    }
    
    // When we get here, either we found a path, or there isn't one.
    return to_return;
}


VCFTraversalFinder::VCFTraversalFinder(const PathHandleGraph& graph, SnarlManager& snarl_manager,
                                       vcflib::VariantCallFile& vcf,
                                       const vector<string>& ref_path_names,
                                       FastaReference* ref_fasta,
                                       FastaReference* ins_fasta,
                                       function<bool(const SnarlTraversal&, int)> skip_alt,
                                       size_t max_traversal_cutoff) :
    graph(graph),
    snarl_manager(snarl_manager),
    skip_alt(skip_alt),
    max_traversal_cutoff(max_traversal_cutoff),
    path_finder(graph, ref_path_names) {

    create_variant_index(vcf, ref_fasta, ins_fasta);
}

VCFTraversalFinder::~VCFTraversalFinder() {
    delete_variant_index();
}

void VCFTraversalFinder::create_variant_index(vcflib::VariantCallFile& vcf, FastaReference* ref_fasta,
                                              FastaReference* ins_fasta) {

    vcflib::Variant var;
#ifdef debug
    cerr << "indexing vcf using alt-path information from graph" << endl;
#endif
    vector<FastaReference*> insertion_fastas;
    if (ins_fasta != nullptr) {
        insertion_fastas.push_back(ins_fasta);
    }
    
    while (vcf.getNextVariant(var)) {
        bool path_found = false;
        path_handle_t path_handle;

        // we need to run this in order for symbolic alleles to get the same hashes as in construct
        if (var.isSymbolicSV()) {
            if (ref_fasta == nullptr) {
                cerr << "[VCFTraversalFinder] Warning: Unable to canonicalize symbolic variant because no reference fasta"
                     << " was given:\n" << var << endl;
                continue;
            }
            bool could_canonicalize = var.canonicalize(*ref_fasta, insertion_fastas, true);
            if (!could_canonicalize) {
                cerr << "[VCFTraversalFinder] Warning: Failed to canonicalize symbolic variant:\n" << var << endl;
                continue;
            }
        }

        // scan paths in the graph for any alt path that could have come from this variant
        // then add any node id from the path to our index
        // we add the first id we find under the assumption that alt paths are entirely contained within sites
        for (int allele = 0; !path_found && allele < var.alleles.size(); ++allele) {
            string alt_path_name = "_alt_" + make_variant_id(var) + "_" + to_string(allele);
            if (graph.has_path(alt_path_name)) {
                path_handle_t path_handle = graph.get_path_handle(alt_path_name);
                if (!graph.is_empty(path_handle)) {
                    path_found = true;
                    step_handle_t step_handle = graph.path_begin(path_handle);
                    handle_t handle = graph.get_handle_of_step(step_handle);
                    id_t node_id = graph.get_id(handle);
                    // copy our variant just this once, and add its new pointer to our map
                    if (node_to_variant.count(node_id)) {
                        node_to_variant[node_id].push_back(new vcflib::Variant(var));
                    } else {
                        node_to_variant[node_id] = list<vcflib::Variant*>({new vcflib::Variant(var)});
                    }
                }
            }
        }
        if (!path_found) {
            cerr << "[VCFTraversalFinder] Warning: No alt path (prefix="
                 << ("_alt_" + make_variant_id(var) + "_") << ") found in graph for variant.  It will be ignored:\n"
                 << var << endl;
        }
    }
#ifdef debug
    cerr << "Indexed " << node_to_variant.size() << " nodes" << endl;
#endif
}

void VCFTraversalFinder::delete_variant_index() {
    for (auto nv : node_to_variant) {
        for (auto var : nv.second) {
            delete var;
        }
    }
    node_to_variant.clear();
}

vector<vcflib::Variant*> VCFTraversalFinder::get_variants_in_site(const Snarl& site) {
    vector<vcflib::Variant*> site_variants;

    pair<unordered_set<id_t>, unordered_set<edge_t> > contents = snarl_manager.deep_contents(&site, graph, false);

    for (auto node_id : contents.first) {
        auto map_it = node_to_variant.find(node_id);
        if (map_it != node_to_variant.end()) {
            for (auto var : map_it->second) {
                site_variants.push_back(var);
            }
        }
    }
            
    return site_variants;
}

pair<vector<pair<SnarlTraversal, vector<int>>>, vector<vcflib::Variant*>>
VCFTraversalFinder::find_allele_traversals(Snarl site) {

    vector<pair<SnarlTraversal, vector<int>>> output_traversals;

    // This traversal finder is pretty simple-minded.  It's expecting forward-oriented variation relative
    // to the reference.  We flip our snarl to canonicalize it if possible.
    if (site.start().backward() && site.end().backward()) {
        Visit start = site.start();
        *site.mutable_start() = site.end();
        *site.mutable_end() = start;
        site.mutable_start()->set_backward(false);
        site.mutable_end()->set_backward(false);
    }
    
    vector<vcflib::Variant*> site_variants = get_variants_in_site(site);

    if (site_variants.empty()) {
        return make_pair(output_traversals, site_variants);
    }

    pair<vector<SnarlTraversal>, vector<pair<step_handle_t, step_handle_t> > > trav_steps = path_finder.find_path_traversals(site);

    if (trav_steps.first.empty()) {
        return make_pair(output_traversals, site_variants);
    }

    // we can certainly relax this if needed
    assert(trav_steps.first.size() == 1);

    step_handle_t start_step = trav_steps.second[0].first;
    step_handle_t end_step = trav_steps.second[0].second;
    path_handle_t ref_path = graph.get_path_handle_of_step(start_step);

#ifdef debug
    cerr << "Computing alt traversals for site " << pb2json(site) << " that spans the following "
         << site_variants.size() << " variants:\n";
    for (auto site_var : site_variants) {
      cerr << " ->" << *site_var << endl;
    }
#endif

    // fill in the alt traversals
    brute_force_alt_traversals(site, site_variants, ref_path, start_step, end_step, output_traversals);

    return make_pair(output_traversals, site_variants);
}

vector<SnarlTraversal> VCFTraversalFinder::find_traversals(const Snarl& site) {
    pair<vector<pair<SnarlTraversal, vector<int>>>, vector<vcflib::Variant*>> allele_travs = find_allele_traversals(site);
    vector<SnarlTraversal> traversals;
    traversals.reserve(allele_travs.first.size());
    for (auto& trav : allele_travs.first) {
        traversals.push_back(trav.first);
    }
    return traversals;
}

void VCFTraversalFinder::brute_force_alt_traversals(
    const Snarl& site,
    const vector<vcflib::Variant*>& site_variants,
    path_handle_t ref_path,
    step_handle_t start_step,
    step_handle_t end_step,
    vector<pair<SnarlTraversal, vector<int> > >& output_traversals) {

    // the haplotype we're going to look for a traversal for
    // it's in terms of alt_alleles below (and not the VCF), so needs
    // to be converted back
    vector<int> haplotype(site_variants.size(), 0);

    // use our skip_alt() method (if defined) to prune the search space
    vector<vector<int>> alt_alleles = get_pruned_alt_alleles(site, site_variants, ref_path);
    assert(alt_alleles.size() == haplotype.size());
    
    // if we failed to prune enough, we print a message here:
    // todo: we can move to a ranking (eg by support), where instead of filtering, we just
    // take the K most supported traversals.  this would avoid ever skipping a site
    if (!check_max_trav_cutoff(alt_alleles)) {
        cerr << "[VCFTraversalFinder] Warning: Site " << pb2json(site) << " with " << site_variants.size()
             << " variants contains too many traversals (>" << max_traversal_cutoff 
             << ") to enumerate so it will be skipped:\n";
        for (auto site_var : site_variants) {
            cerr << "   " << *site_var << endl;
        }
        output_traversals.clear();
        return;
    }

    // increment the haplotype.  we can use this to loop over every possible haplotype
    auto next_haplotype = [&] () -> bool {
        // do this by "adding" 1 to our haplotype. each digit is in base-|alleles|
        for (int i = alt_alleles.size() - 1; i >= 0; --i) {
            if (haplotype[i] < alt_alleles[i].size() - 1) {
                // add to column
                ++haplotype[i];
                return true;
            } else if (i > 0) {
                // carry 1 to left
                haplotype[i] = 0;
            }
        }
        return false;
    };

    do {
        // convert back to vcf allele offsets
        // todo: can we change the enumeration to avoid this?
        vector<int> vcf_haplotype(haplotype.size());
        for (int i = 0; i < site_variants.size(); ++i) {
            vcf_haplotype[i] = (alt_alleles[i][haplotype[i]]);
            assert(skip_alt != nullptr || vcf_haplotype[i] == haplotype[i]);
        }
        // make sure we don't forget the reference.  I'm sure there's a more elegant way to
        // do this, but it's fussy in light of pruning logic
        if (output_traversals.empty() &&
            !std::all_of(vcf_haplotype.begin(), vcf_haplotype.end(), [](int i) {return i == 0;})) {
            vector<int> ref_haplotype(vcf_haplotype.size(), 0);
            pair<SnarlTraversal, bool> alt_traversal = get_alt_traversal(
                site, site_variants, ref_path, start_step, end_step, ref_haplotype);
            assert(alt_traversal.second == true);
            output_traversals.push_back(make_pair(alt_traversal.first, ref_haplotype));
        }
        
        pair<SnarlTraversal, bool> alt_traversal = get_alt_traversal(site, site_variants, ref_path,
                                                                     start_step, end_step, vcf_haplotype);
#ifdef debug
        cerr << "bf haplotype <";
        for (auto allele : vcf_haplotype) {
            cerr << allele << ",";
        }
        cerr << "> gives " << (alt_traversal.second ? "valid" : "invalid") <<  " trav: "
             << pb2json(alt_traversal.first) << endl;
#endif
        if (alt_traversal.second) {
            output_traversals.push_back(make_pair(alt_traversal.first, vcf_haplotype));
        }
    } while (next_haplotype());   
}

pair <SnarlTraversal, bool> VCFTraversalFinder::get_alt_traversal(const Snarl& site,
                                                                  const vector<vcflib::Variant*>& site_variants,
                                                                  path_handle_t ref_path,
                                                                  step_handle_t start_step,
                                                                  step_handle_t end_step,
                                                                  const vector<int>& haplotype) {

    // Find the alt paths that we must cover if we traverse this haplotype
    pair<unordered_set<handle_t>, unordered_set<pair<handle_t, handle_t>>> alt_contents =
        get_haplotype_alt_contents(site_variants, haplotype, ref_path);
    unordered_set<handle_t>& alt_nodes = alt_contents.first;
    unordered_set<pair<handle_t, handle_t>>& alt_edges = alt_contents.second;

    // the edges of our reference path.  we must follow these in our traversal
    // unless we're going to an alt
    unordered_set<pair<handle_t, handle_t> > ref_edges;
    // nodes of the reference path.  we use these to see if we can go from
    // an alt path back to the reference
    unordered_set<handle_t> ref_nodes;
    for (auto step = start_step; step != end_step; step = graph.get_next_step(step)) {
        auto next = graph.get_next_step(step);

        // todo: assuming forward ref path
        ref_edges.insert(graph.edge_handle(graph.get_handle_of_step(step),
                                           graph.get_handle_of_step(next)));

        ref_nodes.insert(graph.get_handle_of_step(step));
    }
    ref_nodes.insert(graph.get_handle_of_step(end_step));

#ifdef debug
    cerr << "  alt nodes: ";
    for (auto alt_node : alt_nodes) {
        cerr << graph.get_id(alt_node) << ":" << graph.get_is_reverse(alt_node) << ",";
    }
    cerr << endl << "  alt edges: ";
    for (auto alt_edge : alt_edges) {
        cerr << graph.get_id(alt_edge.first) << ":" << graph.get_is_reverse(alt_edge.first) << "-"
             << graph.get_id(alt_edge.second) << ":" << graph.get_is_reverse(alt_edge.second) << ",";
    }
    cerr << endl << "  ref edges: ";
    for (auto ref_edge : ref_edges) {
        cerr << graph.get_id(ref_edge.first) << ":" << graph.get_is_reverse(ref_edge.first) << "-"
             << graph.get_id(ref_edge.second) << ":" << graph.get_is_reverse(ref_edge.second) << ",";
    }            
#endif
    
    // we walk by always following reference edges unless we can step into an alt
    // there are some simplifying assumptions about alt paths at play here, like
    // how there are unique hooks between them and the reference.
    bool in_alt_path = false;
    auto walk_forward = [&] (Visit& visit) {
        handle_t handle = graph.get_handle(visit.node_id(), visit.backward());
        // take an alt edge if we find one
        bool found_edge = !graph.follow_edges(handle, false, [&] (const handle_t& next) {
                auto edge = graph.edge_handle(handle, next);
                bool ret = true;
                if (alt_edges.count(edge)) {
                    ret = false; // stop, we found deletion edge
                } else if (alt_nodes.count(next)) {
                    in_alt_path = true;
                    ret = false; // stop, we found an edge to an alt path node
                }
                if (ret == false) {
                    // make sure we never cross this node/edge again in our traversal
                    alt_edges.erase(edge);
                    alt_nodes.erase(graph.get_handle(graph.get_id(next), false));
                    alt_nodes.erase(graph.get_handle(graph.get_id(next), true));
                    visit.set_node_id(graph.get_id(next));
                    visit.set_backward(graph.get_is_reverse(next));
                }
                return ret;
            });
        if (!found_edge) {
            // no alt edge found, take a reference edge
            found_edge = !graph.follow_edges(handle, false, [&] (const handle_t& next) {
                auto edge = graph.edge_handle(handle, next);
                bool ret = true;
                if (ref_edges.count(edge) && ref_nodes.count(next)) {
                    ret = false; // stop, we found a reference edge
                } else if (in_alt_path && ref_nodes.count(next)) {
                    in_alt_path = false;
                    ret = false; // stop, we found a reference node after our alt path
                }
                if (ret == false) {
                    // make sure we never cross this node/edge again in our traversal
                    ref_edges.erase(edge);
                    ref_nodes.erase(graph.get_handle(graph.get_id(next), false));
                    ref_nodes.erase(graph.get_handle(graph.get_id(next), true));
                    visit.set_node_id(graph.get_id(next));
                    visit.set_backward(graph.get_is_reverse(next));
                }
                return ret;
            });
        }
        return found_edge;
    };

    
    Visit visit;
    SnarlTraversal traversal;
    
    // start at the start
    // todo: should make sure this works if our snarl is backward on reference
    visit.set_node_id(graph.get_id(graph.get_handle_of_step(start_step)));
    visit.set_backward(graph.get_is_reverse(graph.get_handle_of_step(start_step)));
    ref_nodes.erase(graph.get_handle(visit.node_id(), false));
    ref_nodes.erase(graph.get_handle(visit.node_id(), true));
    

    if (include_endpoints) {
        *traversal.add_visit() = visit;
    }

#ifdef debug
    cerr << "  start walk: " << pb2json(visit) << endl;
#endif
    
    // walk our traversal
    bool found_end;
    while (walk_forward(visit)) {
#ifdef debug
        cerr << "  visit: " << pb2json(visit) << endl;
#endif
        if (visit.node_id() != graph.get_id(graph.get_handle_of_step(end_step))) {
            *traversal.add_visit() = visit;
        } else {
            found_end = true;
            break;
        }
    }

    if (include_endpoints) {
        Visit* visit = traversal.add_visit();
        handle_t end_handle = graph.get_handle_of_step(end_step);
        visit->set_node_id(graph.get_id(end_handle));
        visit->set_backward(graph.get_is_reverse(end_handle));
    }

    // consistency check: we compare the output to something gotten directly from the
    // path index when doing the reference haplotype.
    if (all_of(haplotype.begin(), haplotype.end(), [] (int i) {return i == 0;})) { 
        SnarlTraversal ref_trav;
        step_handle_t step = graph.get_next_step(start_step);
        if (include_endpoints) {
            Visit* visit = ref_trav.add_visit();
            visit->set_node_id(graph.get_id(graph.get_handle_of_step(start_step)));
            visit->set_backward(graph.get_is_reverse(graph.get_handle_of_step(start_step)));
        }
        for (; step != end_step; step = graph.get_next_step(step)) {
            Visit* visit = ref_trav.add_visit();
            visit->set_node_id(graph.get_id(graph.get_handle_of_step(step)));
            // todo: do we get an orientation out of the path index?  
            visit->set_backward(graph.get_is_reverse(graph.get_handle_of_step(step)));
        }
        if (include_endpoints) {
            Visit* visit = ref_trav.add_visit();
            visit->set_node_id(graph.get_id(graph.get_handle_of_step(end_step)));
            visit->set_backward(graph.get_is_reverse(graph.get_handle_of_step(end_step)));
        }
        assert(found_end && ref_trav == traversal);
    }
    
    return make_pair(traversal, found_end && alt_nodes.empty() && alt_edges.empty());    
}

pair<unordered_set<handle_t>, unordered_set<pair<handle_t, handle_t> > >
VCFTraversalFinder::get_haplotype_alt_contents(
    const vector<vcflib::Variant*>& site_variants,
    const vector<int>& haplotype,
    path_handle_t ref_path) {

    assert(haplotype.size() == site_variants.size());

    unordered_set<handle_t> alt_nodes;
    unordered_set<pair<handle_t, handle_t> > alt_deletion_edges;

    for (size_t allele = 0; allele < haplotype.size(); ++allele) {
        // ignore reference alleles
        if (haplotype[allele] == 0) {
            continue;
        }
        vcflib::Variant* var = site_variants[allele];

        // get the alt path information out of the graph
        pair<SnarlTraversal, vector<edge_t>> alt_path_info = get_alt_path(var, haplotype[allele], ref_path);
        if (alt_path_info.first.visit_size() == 0) {
            // skip deletion alt path where we can't find the deletion edge in the graph
            continue;
        }
        SnarlTraversal& alt_traversal = alt_path_info.first;
        bool is_deletion = !alt_path_info.second.empty();

        if (!is_deletion) {
            // fill in the nodes from the path
            for (size_t i = 0; i < alt_traversal.visit_size(); ++i) {
                alt_nodes.insert(graph.get_handle(alt_traversal.visit(i).node_id(),
                                                  alt_traversal.visit(i).backward()));
            }
        }  else {
            // add the deletion edges from the path
            for (auto deletion_edge : alt_path_info.second) {
                alt_deletion_edges.insert(deletion_edge);
            }
        }
    }

    return make_pair(alt_nodes, alt_deletion_edges);
}

pair<SnarlTraversal, vector<edge_t>> VCFTraversalFinder::get_alt_path(vcflib::Variant* var, int allele,
                                                            path_handle_t ref_path) {

    SnarlTraversal alt_path;
    vector<edge_t> deletion_edges;
    
    string alt_path_name = "_alt_" + make_variant_id(*var) + "_" + to_string(allele);
    if (graph.has_path(alt_path_name) && !graph.is_empty(graph.get_path_handle(alt_path_name))) {
        // if there's an alt path, then we're dealing with a snp or insertion.
        // we take the edges from the path, as well as those back to the reference
        for (handle_t handle : graph.scan_path(graph.get_path_handle(alt_path_name))) {
            // fill in the nodes from the path
            Visit* visit = alt_path.add_visit();
            visit->set_node_id(graph.get_id(handle));
            visit->set_backward(graph.get_is_reverse(handle));
        }
    }  else {
        // there's no alt path, it must be a deletion (if our input allele != 0)
        // in this case we use the reference allele path, and try to find an edge that spans
        // it.  this will be our alt edge
        // todo: put an alt path name maker into utility.hpp
        bool is_deletion = allele != 0;
        alt_path_name = "_alt_" + make_variant_id(*var) + "_0";
        // allele 0 can be empty for an insertion.  we don't complain if it's not in the graph
        assert(allele == 0 || graph.has_path(alt_path_name));

        if (graph.has_path(alt_path_name)) {
            path_handle_t path_handle = graph.get_path_handle(alt_path_name);            
            if (!graph.is_empty(path_handle)) {
                // find where this path begins and ends in the reference path index
                auto first_step_found = step_in_path(graph.get_handle_of_step(graph.path_begin(path_handle)), ref_path);
                assert(first_step_found.second);
                step_handle_t first_step = first_step_found.first;
                auto last_step_found = step_in_path(graph.get_handle_of_step(graph.path_back(path_handle)), ref_path);
                assert(last_step_found.second);
                step_handle_t last_step = last_step_found.first;
                
                // todo: logic needed here if want to support non-forward reference paths.
                first_step = graph.get_previous_step(first_step);
                last_step = graph.get_next_step(last_step);
                if (allele == 0) {
                    handle_t left = graph.get_handle_of_step(first_step);
                    handle_t right = graph.get_handle_of_step(last_step);
                    
                    Visit* visit = alt_path.add_visit();
                    visit->set_node_id(graph.get_id(left));
                    visit->set_backward(graph.get_is_reverse(left));
                    visit = alt_path.add_visit();
                    visit->set_node_id(graph.get_id(right));
                    visit->set_backward(graph.get_is_reverse(right));
                } else {
                    // alt paths don't always line up with deletion edges, so we hunt for
                    // our deletion edge using the path_index here.
                    pair<SnarlTraversal, vector<edge_t>> scan_deletion = scan_for_deletion(var, allele, ref_path,
                                                                                           first_step, last_step);
                    if (scan_deletion.first.visit_size() == 0) {
                        cerr << "[VCFTraversalFinder] Warning: Could not find deletion edge that matches allele "
                             << allele << " of\n" << *var << "\naround alt path" << alt_path_name << ":";
                    }
                    alt_path = std::move(scan_deletion.first);
                    deletion_edges = std::move(scan_deletion.second);
                }
            }
        }
    }

    return make_pair(alt_path, deletion_edges);
}

pair<SnarlTraversal, vector<edge_t>> VCFTraversalFinder::scan_for_deletion(vcflib::Variant* var, int allele, path_handle_t ref_path,
                                                                           step_handle_t first_path_step, step_handle_t last_path_step) {
    assert(allele > 0);

    // if our path matches an edge, we don't need to do anything
    edge_t spanning_edge = graph.edge_handle(graph.get_handle_of_step(first_path_step),
                                             graph.get_handle_of_step(last_path_step));
    if (graph.has_edge(spanning_edge)) {
        SnarlTraversal traversal;
        Visit* visit = traversal.add_visit();
        visit->set_node_id(graph.get_id(graph.get_handle_of_step(first_path_step)));
        visit->set_backward(graph.get_is_reverse(graph.get_handle_of_step(first_path_step)));
        visit = traversal.add_visit();
        visit->set_node_id(graph.get_id(graph.get_handle_of_step(last_path_step)));
        visit->set_backward(graph.get_is_reverse(graph.get_handle_of_step(last_path_step)));
        return make_pair(traversal, vector<edge_t>(1, spanning_edge));
    }

    // we're doing everything via length comparison, so keep track of the length we're
    // looking for (vcf_deletion_length) and the length we have (path_deletion_length)
    int vcf_deletion_length = var->alleles[0].length() - var->alleles[allele].length();

    // make our search window by scanning out the ends of our path
    step_handle_t first_window_step = first_path_step;
    for (int i = 0; i < max_deletion_scan_nodes && graph.has_previous_step(first_window_step); ++i) {
        first_window_step = graph.get_previous_step(first_window_step);
    }
    step_handle_t last_window_step = last_path_step;
    for (int i = 0; i < max_deletion_scan_nodes && graph.has_next_step(last_window_step); ++i) {
        last_window_step = graph.get_next_step(last_window_step);
    }
    
    // index our reference offsets (assuming forward reference path with no cycles)
    // needing this logic doesn't happen very often, otherwise would consider
    // requiring path position interface
    unordered_map<nid_t, int> ref_offsets;
    unordered_map<nid_t, step_handle_t> node_to_step;
    int offset = 0;
    int first_offset = 0;
    int last_offset = 0;
    step_handle_t cur_step = first_window_step;
    while (true) { // want to iterate last-inclusive
        handle_t cur_handle = graph.get_handle_of_step(cur_step);
        assert(graph.get_is_reverse(cur_handle) == false);
        ref_offsets[graph.get_id(cur_handle)] = offset;
        node_to_step[graph.get_id(cur_handle)] = cur_step;
        if (cur_step == first_path_step) {
            first_offset = offset + graph.get_length(cur_handle);
        }
        if (cur_step == last_path_step) {
            last_offset = offset;
        }
        offset += graph.get_length(cur_handle);
        if (cur_step == last_window_step) {
            break;
        } else {
            cur_step = graph.get_next_step(cur_step);
        }
    }

    // find our deletions, and index them by how close they match our given alt path
    // the delta is the min length of the deletion's two endpoints to the paths enpoints. 
    multimap<int, edge_t> delta_to_deletion; // index deleions by their distance from ref path ends
    unordered_map<edge_t, size_t> deletion_to_length; // lookup deletion sizes
    for (cur_step = first_window_step; cur_step != last_window_step; cur_step = graph.get_next_step(cur_step)) {
        handle_t cur_handle = graph.get_handle_of_step(cur_step);
        handle_t next_handle = graph.get_handle_of_step(graph.get_next_step(cur_step));        
        graph.follow_edges(cur_handle, false, [&] (const handle_t& edge_next_handle) {
                if (!graph.get_is_reverse(edge_next_handle) && // ignore inversions
		    graph.get_id(next_handle) != graph.get_id(edge_next_handle) &&
                    ref_offsets.count(graph.get_id(edge_next_handle))) {
                    // we are in a deletion that's contained in the window
                    int deletion_start_offset = ref_offsets[graph.get_id(cur_handle)] + graph.get_length(cur_handle);
                    int deletion_end_offset = ref_offsets[graph.get_id(edge_next_handle)];
                    int delta = std::min(std::abs(deletion_start_offset - first_offset), std::abs(deletion_end_offset - last_offset));
                    delta_to_deletion.insert(make_pair(delta, make_pair(cur_handle, edge_next_handle)));
                    deletion_to_length[make_pair(cur_handle, edge_next_handle)] = deletion_end_offset - deletion_start_offset;
                }
            });
    }

    // our goal is to find a traversal that threads the deletions we find.  in order to do this, our deltions
    // can't overlap
    function<bool(edge_t, const vector<edge_t>&)> doesnt_intersect = [&](edge_t edge, const vector<edge_t>& edge_set) {
        int edge_start = ref_offsets[graph.get_id(edge.first)] + graph.get_length(edge.first);
        int edge_end = ref_offsets[graph.get_id(edge.second)];
        // because of previous assumptins, are edges shoudl always be forward alogn the path.
        // note they are not made wtih graph.edge_handle, and are just oriented in the scan order
        assert(edge_start <= edge_end);
        for (edge_t other_edge : edge_set) {
            int other_start = ref_offsets[graph.get_id(other_edge.first)] + graph.get_length(other_edge.first);
            int other_end = ref_offsets[graph.get_id(other_edge.second)];
            if ((other_start >= edge_start && other_start < edge_end) || (other_end > edge_start && other_end <= edge_end) ||
                (edge_start >= other_start && edge_start < other_end) || (edge_end > other_start && edge_end >= other_end)) {
                return false;
            }
        }
        return true;
    };
    
    // greedily try to find some deletions that add up to the desired length, and are close to spanning
    // the alt path
    int best_delta = numeric_limits<int>::max();
    vector<edge_t> best_set;
    for (auto delta_edge : delta_to_deletion) {
        // can do better than quadratic here, but the sizes should be small enough not to matter
        vector<edge_t> candidate_set = {delta_edge.second};
        size_t total_size = deletion_to_length[delta_edge.second];
        int total_delta = delta_edge.first;
        for (auto delta_edge2 : delta_to_deletion) {
            if (delta_edge2 != delta_edge) {
                if (total_size + deletion_to_length[delta_edge2.second] <= vcf_deletion_length &&
                    // make that cubic...
                    doesnt_intersect(delta_edge2.second, candidate_set)) {
                        total_size += deletion_to_length[delta_edge2.second];
                        total_delta += delta_edge2.first;
                        candidate_set.push_back(delta_edge2.second);
                    }
            }
            if (total_size == vcf_deletion_length) {
                break;
            }
        }
        if (total_delta < best_delta) {
            best_set = candidate_set;
            best_delta = total_delta;
        }
    }

    // sort the edges along the path
    std::sort(best_set.begin(), best_set.end(), [&](edge_t e1, edge_t e2) {
            return ref_offsets[graph.get_id(e1.first)] < ref_offsets[graph.get_id(e2.first)]; });

    SnarlTraversal traversal;
    Visit* visit;
    // fill out the traversal
    for (int i = 0; i < best_set.size(); ++i) {
        // add a visit for each edge endpoint
        if (i == 0 || best_set[i].first != best_set[i-1].second) {
            visit = traversal.add_visit();
            visit->set_node_id(graph.get_id(best_set[i].first));
            visit->set_backward(graph.get_is_reverse(best_set[i].first));
        }
        visit = traversal.add_visit();
        visit->set_node_id(graph.get_id(best_set[i].second));
        visit->set_backward(graph.get_is_reverse(best_set[i].second));
        // the fill in the reference path to the next edge
        if (i < best_set.size() - 1) {
            step_handle_t next_step = node_to_step[graph.get_id(best_set[i + 1].first)];
            step_handle_t cur_step = node_to_step[graph.get_id(best_set[i].second)];
            if (cur_step != next_step) {
                for (cur_step = graph.get_next_step(cur_step); cur_step != next_step; cur_step = graph.get_next_step(cur_step)) {
                    visit = traversal.add_visit();
                    visit->set_node_id(graph.get_id(graph.get_handle_of_step(cur_step)));
                    visit->set_backward(graph.get_is_reverse(graph.get_handle_of_step(cur_step)));
                }
            }
        }
    }

    return make_pair(traversal, best_set);
}


vector<vector<int>> VCFTraversalFinder::get_pruned_alt_alleles(
    const Snarl& site,
    const vector<vcflib::Variant*>& site_variants,
    path_handle_t ref_path) {

    vector<vector<int> > alt_alleles(site_variants.size());

    for (int var_i = 0; var_i < site_variants.size(); ++var_i) {
        for (int allele = 0; allele < site_variants[var_i]->alleles.size(); ++allele) {
            alt_alleles[var_i].push_back(allele);
        }
    }

    // only invoke pruning if we exceed our cutoff.  fairly rare on most graphs
    for (int prune_it = 0; prune_it < max_prune_iterations && !check_max_trav_cutoff(alt_alleles); ++prune_it) {
        for (auto& alleles : alt_alleles) {
            alleles.clear();
        }
        
        for (int var_i = 0; var_i < site_variants.size(); ++var_i) {
            for (int allele = 0; allele < site_variants[var_i]->alleles.size(); ++allele) {
                if (skip_alt == nullptr ||
                    skip_alt(get_alt_path(site_variants[var_i], allele, ref_path).first, prune_it) == false) {
                    alt_alleles[var_i].push_back(allele);
                }
#ifdef debug
                else {
                    cerr << "Pruning allele " << allele << " from variant " << site_variants[var_i]->id << endl;
                }
#endif
            }
            // always leave at least one path through the site, even if that means
            // going through a reference allele that fails the skip_alt check.
            if (alt_alleles[var_i].empty()) {
                alt_alleles[var_i].push_back(0);
            }
        }
    }

    return alt_alleles;
}

bool VCFTraversalFinder::check_max_trav_cutoff(const vector<vector<int> >& alleles) {
    if (alleles.empty()) {
        return true;
    }    
    size_t count = 1;

    for (int i = 0; i < alleles.size(); ++i) {
        count *= alleles[i].size();
        if (count > max_traversal_cutoff) {
            return false;
        }
    }

    return true;
}

pair<step_handle_t, bool> VCFTraversalFinder::step_in_path(handle_t handle, path_handle_t path_handle) const {
    vector<step_handle_t> steps = graph.steps_of_handle(handle);
    // must be a cyclic!
    for (auto step : steps) {
        if (graph.get_path_handle_of_step(step) == path_handle) {
            return make_pair(step, true);
        }
    }
    return make_pair(step_handle_t(), false);
}


FlowTraversalFinder::FlowTraversalFinder(const HandleGraph& graph,
                                         size_t K,
                                         function<double(handle_t)> node_weight_callback,
                                         function<double(edge_t)> edge_weight_callback,
                                         size_t max_traversal_length) :
    graph(graph),
    K(K),
    node_weight_callback(node_weight_callback),
    edge_weight_callback(edge_weight_callback),
    max_traversal_length(max_traversal_length) {
    
}

void FlowTraversalFinder::setK(size_t k) {
    K = k;
}

vector<SnarlTraversal> FlowTraversalFinder::find_traversals(const Snarl& site) {
    return find_weighted_traversals(site).first;
}

vector<Traversal> FlowTraversalFinder::find_traversals(const handle_t& snarl_start, const handle_t& snarl_end) {
    return find_weighted_traversals(snarl_start, snarl_end).first;
}

pair<vector<SnarlTraversal>, vector<double>> FlowTraversalFinder::find_weighted_traversals(const Snarl& site, bool greedy_avg,
                                                                                           const HandleGraph* overlay) {

    const HandleGraph* use_graph = overlay != nullptr ? overlay : & graph;
    
    auto [handle_travs, weights] = find_weighted_traversals(
        use_graph->get_handle(site.start().node_id(), site.start().backward()),
        use_graph->get_handle(site.end().node_id(), site.end().backward()),
        greedy_avg, overlay);

    vector<SnarlTraversal> travs(handle_travs.size());
    for (size_t i = 0; i < handle_travs.size(); ++i) {
        for (const handle_t& h : handle_travs[i]) {
            Visit* visit = travs[i].add_visit();
            visit->set_node_id(use_graph->get_id(h));
            visit->set_backward(use_graph->get_is_reverse(h));
        }
    }

    return make_pair(travs, weights);
}

pair<vector<Traversal>, vector<double>> FlowTraversalFinder::find_weighted_traversals(const handle_t& snarl_start,
                                                                                      const handle_t& snarl_end,
                                                                                      bool greedy_avg,
                                                                                      const HandleGraph* overlay) {

    // Search the overlay graph if one is given.
    const HandleGraph* use_graph = overlay != nullptr ? overlay : & graph;

    size_t max_path_length = max_traversal_length;
    if (max_path_length != numeric_limits<size_t>::max()) {
        // Leave the bounds out of the length limit.
        max_path_length += use_graph->get_length(snarl_start) + use_graph->get_length(snarl_end);
    }
    
    vector<pair<double, vector<handle_t>>> widest_paths = algorithms::yens_k_widest_paths(use_graph, snarl_start, snarl_end, K,
                                                                                          node_weight_callback,
                                                                                          edge_weight_callback,
                                                                                          greedy_avg,
                                                                                          max_path_length);

    vector<Traversal> travs;
    travs.reserve(widest_paths.size());
    vector<double> weights;
    weights.reserve(widest_paths.size());

    for (auto& [weight, path] : widest_paths) {
        weights.push_back(weight);
        travs.emplace_back(std::move(path));
    }

    return make_pair(travs, weights);
}        

GBWTTraversalFinder::GBWTTraversalFinder(const HandleGraph& graph, const gbwt::GBWT& gbwt) : 
    graph(graph),
    gbwt(gbwt) {

}
    
GBWTTraversalFinder::~GBWTTraversalFinder() {

}

pair<vector<SnarlTraversal>, vector<vector<gbwt::size_type>>>
GBWTTraversalFinder::find_gbwt_traversals(const Snarl& site, bool return_paths) {
    
    pair<vector<Traversal>, vector<vector<gbwt::size_type>>> gbwt_travs = this->find_gbwt_traversals(
        graph.get_handle(site.start().node_id(), site.start().backward()),
        graph.get_handle(site.end().node_id(), site.end().backward()),
        return_paths);

    vector<SnarlTraversal> travs(gbwt_travs.first.size());

    for (int64_t i = 0; i < gbwt_travs.first.size(); ++i) {
        SnarlTraversal& proto_trav = travs[i];
        const Traversal& handle_trav = gbwt_travs.first[i];
        for (int64_t j = 0; j < handle_trav.size(); ++j) {
            Visit* visit = proto_trav.add_visit();
            visit->set_node_id(graph.get_id(handle_trav[j]));
            visit->set_backward(graph.get_is_reverse(handle_trav[j]));
        }
    }

    return make_pair(travs, gbwt_travs.second);    
}

pair<vector<Traversal>, vector<vector<gbwt::size_type>>>
GBWTTraversalFinder::find_gbwt_traversals(const handle_t& snarl_start, const handle_t& snarl_end, bool return_paths) {

    // follow all gbwt threads from start to end
    vector<pair<vector<gbwt::node_type>, gbwt::SearchState> > forward_traversals = list_haplotypes(
        graph,
        gbwt,                                                                                                   
        snarl_start,
        [&] (const vector<gbwt::node_type>& new_thread) {
            return gbwt::Node::id(new_thread.back()) == graph.get_id(snarl_end) &&
                gbwt::Node::is_reverse(new_thread.back()) == graph.get_is_reverse(snarl_end);
        });

    // follow all gbwt threads from end to start
    vector<pair<vector<gbwt::node_type>, gbwt::SearchState> > backward_traversals;
    if (!gbwt.bidirectional()) {
        backward_traversals = list_haplotypes(
            graph,
            gbwt,
            graph.flip(snarl_end),
            [&] (const vector<gbwt::node_type>& new_thread) {
                return gbwt::Node::id(new_thread.back()) == graph.get_id(snarl_start) &&
                    gbwt::Node::is_reverse(new_thread.back()) == !graph.get_is_reverse(snarl_start);
            });
    }

    // store them all as snarltraversals
    vector<Traversal> traversals;
    vector<vector<gbwt::size_type>> gbwt_paths;
    traversals.reserve(forward_traversals.size() + backward_traversals.size());

    // copy the forward traversals from gbwt vectors to snarl traversals
    for (int i = 0; i < forward_traversals.size(); ++i) {
        traversals.emplace_back();
        for (auto j = forward_traversals[i].first.begin(); j != forward_traversals[i].first.end(); ++j) {
            traversals.back().push_back(graph.get_handle(gbwt::Node::id(*j), gbwt::Node::is_reverse(*j)));
        }
        if (return_paths) {
            gbwt_paths.push_back(gbwt.locate(forward_traversals[i].second));
        }
    }

    if (!backward_traversals.empty()) {

        // want to check we don't have the same element twice
        std::sort(forward_traversals.begin(), forward_traversals.end(),
                  [&](const pair<vector<gbwt::node_type>, gbwt::SearchState>& t1,
                      const pair<vector<gbwt::node_type>, gbwt::SearchState>& t2) {
                      return t1.first < t2.first; });
        
        // copy and reverse the backward traversals into the snarl traversals
        for (int i = 0; i < backward_traversals.size(); ++i) {

            vector<gbwt::size_type> gbwt_path;
            if (return_paths) {
                gbwt_path = gbwt.locate(backward_traversals[i].second);
            }
            
            // orient along the snarl
            std::reverse(backward_traversals[i].first.begin(), backward_traversals[i].first.end());
            for (auto& gnode : backward_traversals[i].first) {
                gnode = gbwt::Node::encode(gbwt::Node::id(gnode), !gbwt::Node::is_reverse(gnode));
            }

            // search in the forward traversals
            auto si = std::lower_bound(forward_traversals.begin(), forward_traversals.end(), backward_traversals[i],
                                       [&](const pair<vector<gbwt::node_type>, gbwt::SearchState>& t1,
                                           const pair<vector<gbwt::node_type>, gbwt::SearchState>& t2) {
                                           return t1.first < t2.first; });
            if (si != forward_traversals.end() && si->first == backward_traversals[i].first) {
                // we found and exact forward match, just add in the paths
                if (return_paths) {
                    size_t idx = si - forward_traversals.begin();
                    gbwt_paths[idx].insert(gbwt_paths[idx].end(), gbwt_path.begin(), gbwt_path.end());
                }
            } else {
                // insert if not duplicate of existing forward traversal
                traversals.emplace_back();
                for (auto j = backward_traversals[i].first.begin(); j != backward_traversals[i].first.end(); ++j) {
                    traversals.back().push_back(graph.get_handle(gbwt::Node::id(*j), gbwt::Node::is_reverse(*j)));
                }
                if (return_paths) {
                    gbwt_paths.push_back(gbwt.locate(backward_traversals[i].second));
                }
            }
        }
    }
    return make_pair(traversals, gbwt_paths);
}

vector<SnarlTraversal> GBWTTraversalFinder::find_traversals(const Snarl& site) {
    return find_gbwt_traversals(site, false).first;
}

vector<Traversal> GBWTTraversalFinder::find_traversals(const handle_t& snarl_start, const handle_t& snarl_end) {
    return find_gbwt_traversals(snarl_start, snarl_end, false).first;
}

pair<vector<SnarlTraversal>, vector<gbwt::size_type>> GBWTTraversalFinder::find_path_traversals(const Snarl& site) {
    // get the unique traversals
    pair<vector<SnarlTraversal>, vector<vector<gbwt::size_type>>> gbwt_traversals = find_gbwt_traversals(site, true);

    // expand them out to one per path (this is to be consistent with PathTraversalFinder as used in deconstruct)
    pair<vector<SnarlTraversal>, vector<gbwt::size_type>> path_traversals;
    for (size_t i = 0; i < gbwt_traversals.first.size(); ++i) {
        SnarlTraversal& trav = gbwt_traversals.first[i];
        vector<gbwt::size_type>& paths = gbwt_traversals.second[i];
        for (size_t j = 0; j < paths.size(); ++j) {
            path_traversals.first.push_back(trav);
            path_traversals.second.push_back(paths[j]);
        }
    }
    
    return path_traversals;
}

pair<vector<Traversal>, vector<gbwt::size_type>> GBWTTraversalFinder::find_path_traversals(const handle_t& snarl_start, const handle_t& snarl_end) {
    // get the unique traversals
    pair<vector<Traversal>, vector<vector<gbwt::size_type>>> gbwt_traversals = this->find_gbwt_traversals(snarl_start, snarl_end, true);

    // expand them out to one per path (this is to be consistent with PathTraversalFinder as used in deconstruct)
    pair<vector<Traversal>, vector<gbwt::size_type>> path_traversals;
    for (size_t i = 0; i < gbwt_traversals.first.size(); ++i) {
        Traversal& trav = gbwt_traversals.first[i];
        vector<gbwt::size_type>& paths = gbwt_traversals.second[i];
        for (size_t j = 0; j < paths.size(); ++j) {
            path_traversals.first.push_back(trav);
            path_traversals.second.push_back(paths[j]);
        }
    }
    
    return path_traversals;
}



}



