#ifndef VG_GCSA_HELPER_HPP_INCLUDED
#define VG_GCSA_HELPER_HPP_INCLUDED

/** \file 
 * Utility classes and functions for working with GCSA.
 */

#include <gcsa/gcsa.h>
#include <gcsa/lcp.h>

namespace vg {

//------------------------------------------------------------------------------

/*
    These are the proper ways of saving GCSA structures.
    Load them with `vg::io::VPKG::load_one`.
*/

/// Save GCSA to the file.
void save_gcsa(const gcsa::GCSA& index, const std::string& filename, bool show_progress = false);

/// Save LCP array to the file.
void save_lcp(const gcsa::LCPArray& lcp, const std::string& filename, bool show_progress = false);

//------------------------------------------------------------------------------

} // namespace vg

#endif // VG_GCSA_HELPER_HPP_INCLUDED
