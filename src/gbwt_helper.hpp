#ifndef VG_GBWT_HELPER_HPP_INCLUDED
#define VG_GBWT_HELPER_HPP_INCLUDED

/** \file 
 * Utility classes and functions for working with GBWT.
 */

#include <unordered_set>
#include <vector>

#include "position.hpp"

#include <gbwt/dynamic_gbwt.h>
#include <gbwt/fast_locate.h>
#include <handlegraph/mutable_path_handle_graph.hpp>

namespace vg {

std::vector<std::string> parseGenotypes(const std::string& vcf_line, size_t num_samples);

//------------------------------------------------------------------------------

/// Convert gbwt::node_type to handle_t.
inline handle_t gbwt_to_handle(const HandleGraph& graph, gbwt::node_type node) {
    return graph.get_handle(gbwt::Node::id(node), gbwt::Node::is_reverse(node));
}

/// Convert handle_t to gbwt::node_type.
inline gbwt::node_type handle_to_gbwt(const HandleGraph& graph, handle_t handle) {
    return gbwt::Node::encode(graph.get_id(handle), graph.get_is_reverse(handle));
}

/// Convert Mapping to gbwt::node_type.
inline gbwt::node_type mapping_to_gbwt(const Mapping& mapping) {
    return gbwt::Node::encode(mapping.position().node_id(), mapping.position().is_reverse());
}

/// Extract a path as a GBWT path. If the path does not exist, it is treated as empty.
gbwt::vector_type extract_as_gbwt_path(const PathHandleGraph& graph, const std::string& path_name);

/// Find all predecessor nodes of the path, ignoring self-loops. If the path
/// does not exist, it is treated as empty.
gbwt::vector_type path_predecessors(const PathHandleGraph& graph, const std::string& path_name);

//------------------------------------------------------------------------------

// GBWT construction helpers.

/// Determine the node width in bits for the GBWT nodes based on the given graph.
gbwt::size_type gbwt_node_width(const HandleGraph& graph);

//------------------------------------------------------------------------------

/*
    These are the proper ways of saving and loading GBWT structures. In case of
    a failure, these will print an error message and exit the program.

    The vg::io::VPKG interface is effectively the same, but it does not handle
    errors in a consistent way.
*/

/// Load a compressed GBWT from the file.
void load_gbwt(gbwt::GBWT& index, const std::string& filename, bool show_progress = false);

/// Load a dynamic GBWT from the file.
void load_gbwt(gbwt::DynamicGBWT& index, const std::string& filename, bool show_progress = false);

/// Load an r-index from the file.
void load_r_index(gbwt::FastLocate& index, const std::string& filename, bool show_progress = false);

/// Save a compressed GBWT to the file.
void save_gbwt(const gbwt::GBWT& index, const std::string& filename, bool show_progress = false);

/// Save a dynamic GBWT to the file.
void save_gbwt(const gbwt::DynamicGBWT& index, const std::string& filename, bool show_progress = false);

/// Save an r-index to the file.
void save_r_index(const gbwt::FastLocate& index, const std::string& filename, bool show_progress = false);

//------------------------------------------------------------------------------

// FIXME: External GBWT stored in GBZ?
/**
 * Helper class that stores either a GBWT or a DynamicGBWT and loads them from a file
 * or converts between them when necessary.
 */
struct GBWTHandler {
    enum index_type { index_none, index_external, index_compressed, index_dynamic };

    /// Compressed GBWT stored within the handler.
    gbwt::GBWT compressed_owned;

    /// Compressed GBWT stored somewhere else.
    gbwt::GBWT* compressed_external = nullptr;

    /// Dynamic GBWT.
    gbwt::DynamicGBWT dynamic;

    /// Which index is in use.
    index_type in_use = index_none;

    /// The in-memory indexes are backed by this file.
    std::string filename;

    /// Print progress information to stderr when loading/converting indexes.
    bool show_progress = false;

    /// Returns a pointer to the compressed GBWT in use, if any.
    gbwt::GBWT* get_compressed();

    /// Returns a const pointer to the compressed GBWT in use, if any.
    const gbwt::GBWT* get_compressed() const;

    /// Switch to a compressed GBWT, converting it from the dynamic GBWT or reading it
    /// from a file if necessary.
    void use_compressed();

    /// Switch to a dynamic GBWT, converting it from the compressed GBWT or reading it
    /// from a file if necessary.
    void use_dynamic();

    /// Start using this compressed GBWT. Clears the index used as the argument.
    /// The GBWT is not backed by a file.
    void use(gbwt::GBWT& new_index);

    /// Start using this external compressed GBWT. The handler does not take ownership
    /// of the index. This GBWT is not backed by a file.
    void use_external(gbwt::GBWT& new_index);

    /// Start using this dynamic GBWT. Clears the index used as the argument.
    /// The GBWT is not backed by a file.
    void use(gbwt::DynamicGBWT& new_index);

    /// The GBWT is no longer backed by a file.
    void unbacked();

    /// Serialize the in-memory index to this file and start using it as the backing file.
    void serialize(const std::string& new_filename);

    /// Clear the in-memory index.
    void clear();
};

//------------------------------------------------------------------------------

/// Return the list of thread ids / gbwt path ids for the given sample.
std::vector<gbwt::size_type> threads_for_sample(const gbwt::GBWT& gbwt_index, const std::string& sample_name);

/// Insert a GBWT thread into the graph and return its name. Returns an empty string on failure.
/// If a path name is specified and not empty, that name will be used for the inserted path.
/// NOTE: id is a gbwt path id, not a gbwt sequence id.
std::string insert_gbwt_path(MutablePathHandleGraph& graph, const gbwt::GBWT& gbwt_index, gbwt::size_type id, std::string path_name = "");

/// Extract a GBWT thread as a path in the given graph.
/// NOTE: id is a gbwt path id, not a gbwt sequence id.
Path extract_gbwt_path(const HandleGraph& graph, const gbwt::GBWT& gbwt_index, gbwt::size_type id);

//------------------------------------------------------------------------------

/// Copies the reference sample tag from the source GBWT index to the destination GBWT index.
void copy_reference_samples(const gbwt::GBWT& source, gbwt::GBWT& destination);

/// Copies reference samples from the source graph to the destination GBWT index.
/// Every sample with at least one reference path in the source graph is considered
/// a reference sample.
void copy_reference_samples(const PathHandleGraph& source, gbwt::GBWT& destination);

//------------------------------------------------------------------------------

/// Append the nodes of `origin_fragment` onto `last_fragment`, joining the end
/// of a haplotype back to its start. The two vectors must be distinct.
///
/// Linearizing a circular sequence drops a single adjacency, the one connecting
/// the end of the sequence to its start. Appending the origin fragment onto the
/// last fragment restores it. For an unfragmented haplotype the origin and the
/// last fragment are the same path, so it is appended to itself.
void append_wrap_fragment(gbwt::vector_type& last_fragment, const gbwt::vector_type& origin_fragment);

/// Rebuild `source` into a new GBWT in which every haplotype on one of the named
/// `contigs` has its origin (first) fragment appended onto its last fragment
/// (see `append_wrap_fragment`), restoring the adjacency that linearization drops
/// from a circular sequence. All other paths are copied verbatim, and the
/// metadata and the reference sample tag are preserved.
///
/// Fragments are grouped into haplotypes using `gbwt::FragmentMap`, which groups
/// by (sample, contig, phase) and orders fragments by their position along the
/// contig. Every fragment of a haplotype therefore lies on the same contig, and
/// wrapping affects only the last fragment. This means fragmented haplotypes are
/// handled correctly: the wrap joins the true end back to the true start, rather
/// than looping an interior fragment onto itself.
///
/// The origin fragment must contain contig position 0, identified by count 0.
/// VCF-built GBWTs store a dense 0-based fragment identifier in the count field,
/// while GFA W-line GBWTs store the genomic start offset, but in both cases the
/// origin fragment has count 0. If a named contig has a haplotype whose origin
/// fragment is missing (the start was truncated out of the graph), this is an
/// error, as there is no origin node to wrap onto.
///
/// The tail of each wrapped haplotype is assumed to be the true contig end.
/// A truncated end cannot be detected, because the metadata stores no path end
/// coordinate (the GFA W-line end is recomputed as count + length when writing).
/// If the end was truncated, the wrap edge connects the last remaining node to
/// the origin.
gbwt::GBWT wrap_haplotype_paths(const gbwt::GBWT& source, const std::unordered_set<std::string>& contigs);

/// Transform the paths into a GBWT index. Primarily for testing.
gbwt::GBWT get_gbwt(const std::vector<gbwt::vector_type>& paths);

//------------------------------------------------------------------------------

/// Load a translation file (created with vg gbwt --translation) and return a mapping
/// original segment ids to a list of chopped node ids
unordered_map<string, vector<nid_t>> load_translation_map(ifstream& input_stream);

/// Load a translation file (created with vg gbwt --translation) and return a backwards mapping
/// of chopped node to original segment position (id,offset pair)
/// NOTE: hopefully this is just a short-term hack, and we get a general interface baked into
///       the handlegraphs themselves
unordered_map<nid_t, pair<string, size_t>> load_translation_back_map(HandleGraph& graph, ifstream& input_stream);

//------------------------------------------------------------------------------

} // namespace vg

#endif // VG_GBWT_HELPER_HPP_INCLUDED
