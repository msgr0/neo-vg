#include "gcsa_helper.hpp"

#include <iostream>

namespace vg {

//------------------------------------------------------------------------------

void save_gcsa(const gcsa::GCSA& index, const std::string& filename, bool show_progress) {
    if (show_progress) {
        std::cerr << "Saving GCSA to " << filename << std::endl;
    }
    if (!sdsl::store_to_file(index, filename)) {
        std::cerr << "error: [save_gcsa()] cannot write GCSA to " << filename << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

void save_lcp(const gcsa::LCPArray& lcp, const std::string& filename, bool show_progress) {
    if (show_progress) {
        std::cerr << "Saving LCP to " << filename << std::endl;
    }
    if (!sdsl::store_to_file(lcp, filename)) {
        std::cerr << "error: [save_gcsa()] cannot write LCP to " << filename << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

//------------------------------------------------------------------------------

} // namespace vg
