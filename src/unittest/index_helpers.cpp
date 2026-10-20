/** \file
 *
 * Unit tests for utility classes and functions for working with GBWT, GBWTGraph, and GCSA2.
 */

#include "../gbwt_helper.hpp"

#include <gbwtgraph/utils.h>

#include "catch.hpp"



namespace vg {

namespace unittest {

//------------------------------------------------------------------------------

namespace {

// A short haplotype fragment on contig 0 (nodes 11, 12, 14).
gbwt::vector_type wrap_hap_a {
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(11, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(12, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(14, false)),
};

// A short haplotype fragment on contig 1 (nodes 21, 22, 24).
gbwt::vector_type wrap_hap_b {
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(21, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(22, false)),
    static_cast<gbwt::vector_type::value_type>(gbwt::Node::encode(24, false)),
};

// Extract the forward path for a given path id from a bidirectional GBWT.
gbwt::vector_type extract_forward(const gbwt::GBWT& index, gbwt::size_type path_id) {
    return index.extract(gbwt::Path::encode(path_id, false));
}

} // anonymous namespace

TEST_CASE("append_wrap_fragment appends the origin onto the last fragment", "[index_helpers]") {
    SECTION("distinct fragments") {
        gbwt::vector_type last = wrap_hap_b;
        append_wrap_fragment(last, wrap_hap_a);
        gbwt::vector_type expected = wrap_hap_b;
        expected.insert(expected.end(), wrap_hap_a.begin(), wrap_hap_a.end());
        REQUIRE(last == expected);
    }

    SECTION("unfragmented haplotype is appended to itself") {
        gbwt::vector_type path = wrap_hap_a;
        gbwt::vector_type origin = wrap_hap_a;
        append_wrap_fragment(path, origin);
        gbwt::vector_type expected = wrap_hap_a;
        expected.insert(expected.end(), wrap_hap_a.begin(), wrap_hap_a.end());
        REQUIRE(path == expected);
    }
}

TEST_CASE("wrap_haplotype_paths wraps unfragmented haplotypes", "[index_helpers]") {
    // Two contigs, each with a single haplotype (sample "sampleA"), plus a
    // generic reference path on contig 0 that must never be wrapped.
    std::vector<gbwt::vector_type> source {
        wrap_hap_a, // path 0: haplotype on contig 0, count 0
        wrap_hap_b, // path 1: haplotype on contig 1, count 0
        wrap_hap_a, // path 2: generic reference on contig 0, count 0
    };
    gbwt::GBWT index = get_gbwt(source);
    index.addMetadata();
    index.metadata.setSamples(std::vector<std::string>{ "sampleA", gbwtgraph::GENERIC_PATH_SAMPLE_NAME });
    index.metadata.setContigs(std::vector<std::string>{ "chr1", "chrM" });
    index.metadata.setHaplotypes(1);
    index.metadata.addPath(0, 0, 0, 0); // sampleA, chr1, phase 0, count 0
    index.metadata.addPath(0, 1, 0, 0); // sampleA, chrM, phase 0, count 0
    index.metadata.addPath(1, 0, 0, 0); // generic,  chr1, phase 0, count 0

    SECTION("wrap one contig") {
        std::unordered_set<std::string> contigs { "chrM" };
        gbwt::GBWT wrapped = wrap_haplotype_paths(index, contigs);

        // chr1 haplotype untouched.
        REQUIRE(extract_forward(wrapped, 0) == wrap_hap_a);
        // The chrM haplotype is unfragmented, so wrapping appends it to itself.
        gbwt::vector_type expected_b = wrap_hap_b;
        expected_b.insert(expected_b.end(), wrap_hap_b.begin(), wrap_hap_b.end());
        REQUIRE(extract_forward(wrapped, 1) == expected_b);
        // Generic reference on chr1 untouched.
        REQUIRE(extract_forward(wrapped, 2) == wrap_hap_a);
        // Metadata preserved.
        REQUIRE(wrapped.metadata == index.metadata);
    }

    SECTION("wrapping a generic reference contig does not affect the reference") {
        // chr1 has both a haplotype and a generic path; wrapping chr1 should
        // wrap only the haplotype, not the generic reference path.
        std::unordered_set<std::string> contigs { "chr1" };
        gbwt::GBWT wrapped = wrap_haplotype_paths(index, contigs);

        gbwt::vector_type expected_a = wrap_hap_a;
        expected_a.insert(expected_a.end(), wrap_hap_a.begin(), wrap_hap_a.end());
        REQUIRE(extract_forward(wrapped, 0) == expected_a); // haplotype wrapped
        REQUIRE(extract_forward(wrapped, 1) == wrap_hap_b); // chrM untouched
        REQUIRE(extract_forward(wrapped, 2) == wrap_hap_a); // generic untouched
    }

    SECTION("empty contig set is a no-op copy") {
        std::unordered_set<std::string> contigs;
        gbwt::GBWT wrapped = wrap_haplotype_paths(index, contigs);
        REQUIRE(extract_forward(wrapped, 0) == wrap_hap_a);
        REQUIRE(extract_forward(wrapped, 1) == wrap_hap_b);
        REQUIRE(extract_forward(wrapped, 2) == wrap_hap_a);
    }
}

TEST_CASE("wrap_haplotype_paths wraps fragmented haplotypes onto the last fragment", "[index_helpers]") {
    // A haplotype on chrM split into two fragments: an origin fragment
    // (count 0) and a second fragment (count 1). Wrapping must append the
    // origin onto the last fragment, leaving the origin fragment itself
    // untouched, so that the single wrap adjacency joins the true end back to
    // the true start.
    std::vector<gbwt::vector_type> source {
        wrap_hap_a, // path 0: origin fragment on chrM, count 0
        wrap_hap_b, // path 1: last fragment on chrM,   count 1
    };
    gbwt::GBWT index = get_gbwt(source);
    index.addMetadata();
    index.metadata.setSamples(std::vector<std::string>{ "sampleA" });
    index.metadata.setContigs(std::vector<std::string>{ "chrM" });
    index.metadata.setHaplotypes(1);
    index.metadata.addPath(0, 0, 0, 0); // sampleA, chrM, phase 0, count 0 (origin)
    index.metadata.addPath(0, 0, 0, 1); // sampleA, chrM, phase 0, count 1 (last)

    std::unordered_set<std::string> contigs { "chrM" };
    gbwt::GBWT wrapped = wrap_haplotype_paths(index, contigs);

    // Origin fragment is copied verbatim, not doubled.
    REQUIRE(extract_forward(wrapped, 0) == wrap_hap_a);
    // Last fragment has the origin appended onto it.
    gbwt::vector_type expected_last = wrap_hap_b;
    expected_last.insert(expected_last.end(), wrap_hap_a.begin(), wrap_hap_a.end());
    REQUIRE(extract_forward(wrapped, 1) == expected_last);
}

TEST_CASE("wrap_haplotype_paths preserves reference samples", "[index_helpers]") {
    std::vector<gbwt::vector_type> source {
        wrap_hap_a, // path 0: reference sample on contig 0
        wrap_hap_b, // path 1: haplotype on contig 1
    };
    gbwt::GBWT index = get_gbwt(source);
    index.addMetadata();
    index.metadata.setSamples(std::vector<std::string>{ "refSample", "sampleA" });
    index.metadata.setContigs(std::vector<std::string>{ "chr1", "chrM" });
    index.metadata.setHaplotypes(2);
    index.metadata.addPath(0, 0, 0, 0); // refSample, chr1
    index.metadata.addPath(1, 1, 0, 0); // sampleA,   chrM
    index.tags.set(gbwtgraph::REFERENCE_SAMPLE_LIST_GBWT_TAG, "refSample");

    std::unordered_set<std::string> contigs { "chr1", "chrM" };
    gbwt::GBWT wrapped = wrap_haplotype_paths(index, contigs);

    // The reference-sample path on chr1 has REFERENCE sense, so it must not be
    // wrapped even though chr1 is named.
    REQUIRE(extract_forward(wrapped, 0) == wrap_hap_a);
    // The unfragmented haplotype on chrM is appended to itself.
    gbwt::vector_type expected_b = wrap_hap_b;
    expected_b.insert(expected_b.end(), wrap_hap_b.begin(), wrap_hap_b.end());
    REQUIRE(extract_forward(wrapped, 1) == expected_b);
    // Reference-samples tag preserved.
    REQUIRE(wrapped.tags.get(gbwtgraph::REFERENCE_SAMPLE_LIST_GBWT_TAG) == "refSample");
}

TEST_CASE("wrap_haplotype_paths rejects a missing origin fragment", "[index_helpers]") {
    // A single haplotype fragment on contig chrM whose count is 1 (its origin,
    // count 0, was truncated out of the graph). Wrapping must throw.
    std::vector<gbwt::vector_type> source { wrap_hap_b };
    gbwt::GBWT index = get_gbwt(source);
    index.addMetadata();
    index.metadata.setSamples(std::vector<std::string>{ "sampleA" });
    index.metadata.setContigs(std::vector<std::string>{ "chrM" });
    index.metadata.setHaplotypes(1);
    index.metadata.addPath(0, 0, 0, 1); // count 1: no origin fragment

    std::unordered_set<std::string> contigs { "chrM" };
    REQUIRE_THROWS_AS(wrap_haplotype_paths(index, contigs), std::runtime_error);
}

//------------------------------------------------------------------------------

}
}
