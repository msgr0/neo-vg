/// \file append_path_move.cpp
///
/// Unit tests for the append_path() overload that takes its second path by
/// rvalue reference.

#include "../../src/unittest/catch.hpp"

#include "../../src/path.hpp"
#include "vg/io/json2pb.h"

namespace vg {
namespace unittest {

TEST_CASE("append_path takes the mappings of a temporary path", "[path][append_path]") {
    // Build a path that visits the given nodes with one 4 bp match each.
    auto make_path = [](const std::vector<int64_t>& node_ids) {
        Path path;
        for (int64_t node_id : node_ids) {
            Mapping* mapping = path.add_mapping();
            mapping->mutable_position()->set_node_id(node_id);
            Edit* edit = mapping->add_edit();
            edit->set_from_length(4);
            edit->set_to_length(4);
        }
        return path;
    };

    Path copied_onto = make_path({1});
    Path moved_onto = make_path({1});
    Path source = make_path({2, 3});

    append_path(copied_onto, source);
    append_path(moved_onto, std::move(source));

    SECTION("the result matches appending a copy") {
        REQUIRE(moved_onto.mapping_size() == 3);
        REQUIRE(pb2json(moved_onto) == pb2json(copied_onto));
    }

    SECTION("the source path has no mappings left") {
        REQUIRE(source.mapping_size() == 0);
    }
}

}
}
