#ifndef VG_TRAVERSAL_FINDER_HPP_INCLUDED
#define VG_TRAVERSAL_FINDER_HPP_INCLUDED
// traversal_finder.hpp: TraversalFinder interfaces

// TraversalFinder is an interface for finding paths through Snarls.  The various
// subclasses implement different algorithms for doing this.  A traversal of a snarl
// can be thought of as an allele in a site.

#include <iostream>
#include <algorithm>
#include <functional>
#include <cmath>
#include <regex>
#include <limits>
#include <unordered_set>
#include <unordered_map>
#include <list>
#include <limits>

#include <structures/immutable_list.hpp>

#include <vg/vg.pb.h>
#include "vg.hpp"
#include "hash_map.hpp"
#include "utility.hpp"
#include "types.hpp"
#include "snarls.hpp"
#include "path_index.hpp"
#include "genotypekit.hpp"
#include "gbwt_helper.hpp"

namespace vg {

using namespace std;


// some Protobuf replacements
using Traversal = vector<handle_t>;
using PathInterval = pair<step_handle_t, step_handle_t>;
string traversal_to_string(const PathHandleGraph* graph, const Traversal& traversal, int64_t max_steps = 10);
// replaces pb2json(snarl)
string graph_interval_to_string(const HandleGraph* graph, const handle_t& start_handle, const handle_t& end_handle);
string path_interval_to_string(const PathHandleGraph* graph, const PathInterval& path_interval);

/// A finder that lists walks through a site, given the handles of its bounds.
class SiteTraversalFinder {
public:
    virtual ~SiteTraversalFinder() = default;

    /**
     * Find walks through the site between two bound handles.
     *
     * @param snarl_start The start bound of the site, oriented into the site.
     * @param snarl_end The end bound of the site, oriented out of the site.
     * @return Walks that each start with snarl_start and end with snarl_end.
     */
    virtual vector<Traversal> find_traversals(const handle_t& snarl_start, const handle_t& snarl_end) = 0;
};

/**
 * Represents a strategy for finding traversals of (nested) sites. Polymorphic
 * base class/interface.
 */
class TraversalFinder {
public:
    virtual ~TraversalFinder() = default;
    
    virtual vector<SnarlTraversal> find_traversals(const Snarl& site) = 0;
};

class ExhaustiveTraversalFinder : public TraversalFinder {
    
    const HandleGraph& graph;
    SnarlManager& snarl_manager;
    bool include_reversing_traversals;
    
public:
    ExhaustiveTraversalFinder(const HandleGraph& graph, SnarlManager& snarl_manager,
                              bool include_reversing_traversals = false);
    
    virtual ~ExhaustiveTraversalFinder();
    
    /**
     * Exhaustively enumerate all traversals through the site. Only valid for
     * acyclic Snarls.
     */
    virtual vector<SnarlTraversal> find_traversals(const Snarl& site);
    
protected:
    void stack_up_valid_walks(handle_t walk_head, vector<Visit>& stack);
    virtual bool visit_next_node(handle_t handle) { return true; }
    void add_traversals(vector<SnarlTraversal>& traversals, handle_t traversal_start,
                        unordered_set<handle_t>& stop_at, unordered_set<handle_t>& yield_at);
};


/** This finder
 * uses the PathHandleGraph interface.
 * It returns all traversals through a snarl that are contained
 * within paths in the graph.  It can also return a mapping from the traversals
 * to their paths*/
class PathTraversalFinder : public TraversalFinder, public SiteTraversalFinder {
    
protected:
    // our graph with indexed path positions
    const PathHandleGraph& graph;
    
    // restrict to these paths
    unordered_set<path_handle_t> paths;
    
public:
    // if path_names not empty, only those paths will be considered
    PathTraversalFinder(const PathHandleGraph& graph,
                        const vector<string>& path_names = {});

    /**
     * Return all traversals through the site that are sub-paths of embedded paths in the graph
     */
    virtual vector<SnarlTraversal> find_traversals(const Snarl& site);
    vector<Traversal> find_traversals(const handle_t& snarl_start, const handle_t& snarl_end) override;

    /**
    * Like above, but return the path steps for the for the traversal endpoints
    */
    virtual pair<vector<SnarlTraversal>, vector<PathInterval>> find_path_traversals(const Snarl& site);
    virtual pair<vector<Traversal>, vector<PathInterval>> find_path_traversals(const handle_t& snarl_start, const handle_t& snarl_end);

};    
    

/**
 * This traversal finder finds one or more traversals through leaf sites with
 * no children. It uses a depth-first search. It doesn't work on non-leaf
 * sites, and is not guaranteed to find all traversals. Only works on acyclic
 * sites that are start-end-reachable.
 */
class TrivialTraversalFinder : public TraversalFinder {

    // Holds the vg graph we are looking for traversals in.
    const HandleGraph& graph;

public:
    TrivialTraversalFinder(const HandleGraph& graph);

    virtual ~TrivialTraversalFinder() = default;
    
    /**
     * Find at least one traversal of the site by depth first search, if any
     * exist. Only works on sites with no children.
     */
    virtual vector<SnarlTraversal> find_traversals(const Snarl& site);
};

/**
 * This TraversalFinder returns a traversals and their corresponding genotypes
 * from an input vcf. It relies on alt-paths in the graph (via construct -a)
 * to map between the vcf and the graph.  
 */
class VCFTraversalFinder : public TraversalFinder {

protected:
    const PathHandleGraph& graph;

    /// Use this to check if our snarl runs through a reference path
    /// (may be overkill, but can be used for debugging)
    PathTraversalFinder path_finder;
    
    /// The SnarlManager managiung the snarls we use
    SnarlManager& snarl_manager;

    /// Store variants indexed by an arbitrary node in one of their associated
    /// alt paths. We can then use this to find all variants in a top-level snarl
    unordered_map<id_t, list<vcflib::Variant*>> node_to_variant;

    /// Use this method to prune the search space by selecting alt-alleles
    /// to skip by considering their paths (in SnarlTraversal) format
    /// It will try again and again until enough traversals are pruned,
    /// with iteration keeping track of how many tries (so it should become stricter
    /// as iteration increases)
    function<bool(const SnarlTraversal& alt_path, int iteration)> skip_alt;

    /// If a snarl has more than this many traversals, return nothing and print
    /// a warning.  Dense and large deletions will make this happen from time
    /// to time.  In practice, skip_alt (above) can be used to prune down
    /// the search space by selecting alleles to ignore.
    size_t max_traversal_cutoff;

    /// Maximum number of pruning iterations
    size_t max_prune_iterations = 2;
    
    /// Include snarl endpoints in traversals
    bool include_endpoints = true;

    /// How far to scan when looking for deletions
    size_t max_deletion_scan_nodes = 50;

public:

    /**
     * Make a new VCFTraversalFinder.  Builds the indexes needed to find all the 
     * variants in a site.
     *
     * The skip_alt() method is defined, it is run on the alt-path of each variant
     * allele in the snarl.  If it returns true, that alt-path will never be included
     * in any traversals returned in find_traversals().  
     * This is used to, for example, use read support to prune the number of traversals
     * that are enumerated.
     */
    VCFTraversalFinder(const PathHandleGraph& graph, SnarlManager& snarl_manager, vcflib::VariantCallFile& vcf,
                       const vector<string>& ref_path_names = {},
                       FastaReference* fasta_ref = nullptr,
                       FastaReference* ins_ref = nullptr,
                       function<bool(const SnarlTraversal&, int)> skip_alt = nullptr,
                       size_t max_traversal_cutoff = 50000);
        
    virtual ~VCFTraversalFinder();
    
    /**
     * Find traversals for the site.  Each traversa is returned in a pair with
     * its haplotype.  The haplotype refers to the list of variants (also returned)
     */
    pair<vector<pair<SnarlTraversal, vector<int>>>, vector<vcflib::Variant*>> find_allele_traversals(Snarl site);

    /**
     * Return a list of traversals for the site.  The same traversals as above, only the
     * haplotype information not included
     */
    virtual vector<SnarlTraversal> find_traversals(const Snarl& site);

    /** 
     * Get all the variants that are contained in a site */
    vector<vcflib::Variant*> get_variants_in_site(const Snarl& site);

    
protected:

    /** Load up all the variants into our node index
     */
    void create_variant_index(vcflib::VariantCallFile& vcf, FastaReference* ref_fasta = nullptr,
                              FastaReference* ins_fasta = nullptr);
    void delete_variant_index();

    /** Get a traversal for every possible haplotype (but reference)
     * in the most naive way possibe.  This will blow up terribly for sites that contain more than a few
     * variants.  There's an obvious dynamic programming speedup, but the main issue is that
     * the output size is exponential in the number of variants.
     */
    void brute_force_alt_traversals(const Snarl& site,
                                    const vector<vcflib::Variant*>& site_variants,
                                    path_handle_t ref_path,
                                    step_handle_t start_step,
                                    step_handle_t end_step,
                                    vector<pair<SnarlTraversal, vector<int> > >& output_traversals);

    /** Get a traversal for a given haplotype.  It gets all the nodes and edges from the alt 
     * paths, and greedily walks over them whenever possible (traversing the reference otherwise).
     * if there is no traversal that can satisfy the haplotype, then the returned bool is set to false
     */
    pair<SnarlTraversal, bool> get_alt_traversal(const Snarl& site,
                                                 const vector<vcflib::Variant*>& site_variants,
                                                 path_handle_t ref_path,
                                                 step_handle_t start_step,
                                                 step_handle_t end_step,
                                                 const vector<int>& haplotype);

    /** Get a set of all alt path nodes and deletion edges for a halptype.
     */
    pair<unordered_set<handle_t>, unordered_set<pair<handle_t, handle_t> >>
         get_haplotype_alt_contents(const vector<vcflib::Variant*>& site_variants,
                                    const vector<int>& haplotype,
                                    path_handle_t ref_path);

    /** Get one alt-path out of the graph in the form of a snarl traversal.  if the path is a deletion,
     *  the edges corresponding to the deletion are also returned.  note that it is indeed possible
     *  for one alt path (and therefore one vcf alleles) to correspond to several deletion edges in the 
     *  graph due to normalization during construction.
     */
    pair<SnarlTraversal, vector<edge_t>> get_alt_path(vcflib::Variant* site_variant, int allele, path_handle_t ref_path);

    /**
     * An alt path for a deletion is the deleted reference path.  But sometimes vg construct doesn't
     * write a deletion edge that exactly jumps over the alt path.  In these cases, we need to 
     * search the graph for one. This does a brute-force check of all deletion edges in the vicinity
     * for one that's the same size as the one we're looking for.  
     * It tries to find a set of nearyby deletions that match the desired length.
     * Todo: check the sequence as well
     * Also todo: It'd be really nice if construct -fa would make the deletion-edge easily inferrable 
     * from the alt path.  It really shouldn't be necessary to hunt around. 
     * Returns: <deletion traversal, list of deletion edges>
     */
    pair<SnarlTraversal, vector<edge_t>> scan_for_deletion(vcflib::Variant* var, int allele, path_handle_t ref_path,
                                                           step_handle_t first_path_step, step_handle_t last_path_step);

    /**
     * Prune our search space using the skip_alt method.  Will return a list of pruned VCF alleles/
     *
     * ex, if the input has A --> T
     *                      G --> C,A
     * there input alleles are <0,1>, <0,1,2>.  If there's no support for the G->C on the second one,
     * the output would be <0,1>, <0,2>.
     *
     */
    vector<vector<int>> get_pruned_alt_alleles(const Snarl& site,
                                               const vector<vcflib::Variant*>& site_variants,
                                               path_handle_t ref_path);

    /**
     * Count the possible traversal paths.  Return false if we ever get beyond our cutoff
     */
    bool check_max_trav_cutoff(const vector<vector<int> >& alleles);

    /**
     * Lookup a node in the reference path (mimics old PathIndex)
     */
    pair<step_handle_t, bool> step_in_path(handle_t handle, path_handle_t path_handle) const;
                                
};

/** Finds traversals with the most flow.  Node and edge weights are specified 
 * using the callbacks and can be used, ex, to yield read supports.  
 * If one traversal is requested, then the path with the highest flow (whose
 * node or edge with the minimum weight is maximum) is returned.  If K
 * traversals are specified, then the K highest flow traversals are returned.
 * This is designed to be a replacement for RepresentativeTraversalFinder.
 * It should do a better job of enumerating off-reference traversals, and will
 * of course guarantee to return all the optimal traversals (in the context of max flow).
 * Unlike RepresentativeTraversalFinder, it does not currently support nested
 * snarls, so all traversals returned are explicit.  
 * It is possible that it will blow up on massive snarls, espeically for large Ks. 
 */ 
class FlowTraversalFinder : public TraversalFinder, public SiteTraversalFinder {
    
protected:
    /// The graph to search, unless a search is given an overlay graph to use
    /// instead.
    const HandleGraph& graph;

    /// The K-best traversals are returned
    size_t K;

    /// Callbacks to get supports
    function<double(handle_t)> node_weight_callback;
    function<double(edge_t)> edge_weight_callback;

    /// Call off the search as soon as a traversal of this length (bp) is encountered
    size_t max_traversal_length;
    
public:
    
    /**
     * Make a finder that searches a graph for the widest traversals.
     *
     * @param graph The graph to search. The finder keeps a reference to it,
     * so it must outlive the finder.
     * @param K The most traversals to find for each site.
     * @param node_weight_callback A function giving the flow through a node.
     * @param edge_weight_callback A function giving the flow through an edge.
     * @param max_traversal_length A length in bases. When a search finds a
     * traversal longer than this, not counting the site's bounds, it stops
     * and returns the traversals found so far, including that one.
     */
    FlowTraversalFinder(const HandleGraph& graph,
                        size_t K,
                        function<double(handle_t)> node_weight_callback,
                        function<double(edge_t)> edge_weight_callback,
                        size_t max_traversal_length = numeric_limits<size_t>::max());

    /**
     * Find the K widest (most flow) traversals through a site.
     *
     * @return The traversals, widest first.
     */
    virtual vector<SnarlTraversal> find_traversals(const Snarl& site);

    /**
     * Find the K widest traversals through the site between two bound
     * handles.
     *
     * @return The traversals, widest first.
     */
    vector<Traversal> find_traversals(const handle_t& snarl_start, const handle_t& snarl_end) override;

    /**
     * Find the K widest traversals through a site, with their flows.
     *
     * @param greedy_avg True to search heuristically by the running average
     * flow, which may miss the widest traversals, rather than by the
     * smallest flow along each traversal.
     * @param overlay A graph to search instead of the finder's graph, or null.
     * @return The traversals, and a vector whose i-th entry is the flow
     * through the i-th traversal. Without greedy_avg, the traversals come
     * widest first. With greedy_avg, each flow is an average of the weights
     * along the traversal, and a traversal can have more flow than the one
     * before it.
     */
    virtual pair<vector<SnarlTraversal>, vector<double>> find_weighted_traversals(const Snarl& site,
                                                                                  bool greedy_avg = false,
                                                                                  const HandleGraph* overlay = nullptr);

    /**
     * Find the K widest traversals through the site between two bound
     * handles, with their flows.
     *
     * @param snarl_start The start bound of the site, oriented into the site.
     * @param snarl_end The end bound of the site, oriented out of the site.
     * @param greedy_avg True to search heuristically by the running average
     * flow, which may miss the widest traversals, rather than by the
     * smallest flow along each traversal.
     * @param overlay A graph to search instead of the finder's graph, or null.
     * If given, the bound handles must be in the overlay, and the traversals
     * are returned as handles in it.
     * @return The traversals, and a vector whose i-th entry is the flow through the i-th traversal. Without
     * greedy_avg, the traversals come widest first. With greedy_avg, each
     * flow is an average of the weights along the traversal, and a traversal
     * can have more flow than the one before it.
     */
    virtual pair<vector<Traversal>, vector<double>> find_weighted_traversals(const handle_t& snarl_start,
                                                                             const handle_t& snarl_end,
                                                                             bool greedy_avg = false,
                                                                             const HandleGraph* overlay = nullptr);

    /// Set K
    void setK(size_t k);

};    

/** Rerturn all traversals of a snarl that correspond to haplotypes stored in a GBWT
 */
class GBWTTraversalFinder : public TraversalFinder, public SiteTraversalFinder {

protected:
    
    const HandleGraph& graph;
    const gbwt::GBWT& gbwt;
    
public:
    
    GBWTTraversalFinder(const HandleGraph& graph, const gbwt::GBWT& gbwt);
    
    virtual ~GBWTTraversalFinder();

    /* Return a traversal for every gbwt thread through the snarl 
     */
    virtual vector<SnarlTraversal> find_traversals(const Snarl& site);
    vector<Traversal> find_traversals(const handle_t& snarl_start, const handle_t& snarl_end) override;

    /** Return the traversals, paired with their path identifiers in the gbwt.  The traversals are 
     *  unique, but there can be more than one path along each one (hence the vector)
     */
    virtual pair<vector<SnarlTraversal>, vector<vector<gbwt::size_type>>>
    find_gbwt_traversals(const Snarl& site, bool return_paths = true);
    virtual pair<vector<Traversal>, vector<vector<gbwt::size_type>>>
    find_gbwt_traversals(const handle_t& snarl_start, const handle_t& snarl_end, bool return_paths = true);

    /** Return traversals paired with path identifiers from the GBWT.  The traversals are *not* unique
     * (which is consistent with PathTraversalFinder)
     * To get the sample name from the path identifier id, use gbwtgraph::get_path_sample_name();
     */
    virtual pair<vector<SnarlTraversal>, vector<gbwt::size_type>> find_path_traversals(const Snarl& site);
    virtual pair<vector<Traversal>, vector<gbwt::size_type>> find_path_traversals(const handle_t& snarl_start, const handle_t& snarl_end);

    const gbwt::GBWT& get_gbwt() { return gbwt; }
    
protected:

    /**
     * Breadth first search from the start to the end, only branching if there's a haplotype 
     * in the GBWT, and returning all unique haplotypes found. 
     */
    vector<pair<vector<gbwt::node_type>, gbwt::SearchState> > get_spanning_haplotypes(handle_t start, handle_t end);    
};

}

#endif
