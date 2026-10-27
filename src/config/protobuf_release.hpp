#ifndef VG_CONFIG_PROTOBUF_RELEASE_HPP_INCLUDED
#define VG_CONFIG_PROTOBUF_RELEASE_HPP_INCLUDED

// Force-included into every vg and libvgio C++ file (see the Makefile).
//
// We build without NDEBUG so that assert() stays on. But protobuf 6 keys its
// debug hardening on NDEBUG too: a has-bit consistency check in every message
// destructor, merge and serialize, and forced copies of default strings. That
// cost 11% of instructions on HiFi reads. So include protobuf's headers once
// with NDEBUG set (their include guards then keep the release behavior), and
// put assert() back on afterwards.

// Files built without protobuf on the include path (like CMake's compiler check) skip this.
#if !defined(NDEBUG) && __has_include(<google/protobuf/port.h>)
#define NDEBUG
#include <google/protobuf/port.h>
#include <google/protobuf/repeated_field.h>
#include <google/protobuf/repeated_ptr_field.h>
#include <google/protobuf/map.h>
#include <google/protobuf/generated_message_tctable_impl.h>
#undef NDEBUG
// assert.h has no include guard for assert itself; re-including it redefines assert() as enabled.
#include <assert.h>
#endif

#endif
