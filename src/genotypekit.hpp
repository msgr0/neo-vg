#ifndef VG_GENOTYPEKIT_HPP_INCLUDED
#define VG_GENOTYPEKIT_HPP_INCLUDED
// genotypekit.hpp: defines pluggable modules for building the genotyper

// The basic idea here is we're going to create a few of these classes, fill in
// their public parameter fields, and then wire them up and set them going to
// emit genotypes in a streaming fashion.

#include <iostream>
#include <algorithm>
#include <functional>
#include <cmath>
#include <regex>
#include <limits>
#include <unordered_set>
#include <unordered_map>
#include <list>
#include <vg/vg.pb.h>
#include "vg.hpp"
#include "hash_map.hpp"
#include "types.hpp"
#include "statistics.hpp"
#include "snarls.hpp"
#include "path_index.hpp"
#include "packer.hpp"

namespace vg {

using namespace std;

// We also supply utility functions for working with genotyping Protobuf objects

/**
 * Get the total read support in a Support.
 */
double total(const Support& support);

/**
 * Get the minimum support of a pair of Supports, by taking the min in each
 * orientation.
 */
Support support_min(const Support& a, const Support& b);

/**
 * Get the maximum support of a pair of Supports, by taking the max in each
 * orientation.
 */
Support support_max(const Support& a, const Support& b);

/**
 * Flip the orientations of a Support.
 */
Support flip(const Support& to_flip);

/**
 * Add two Support values together, accounting for strand.
 */
Support operator+(const Support& one, const Support& other);

/**
 * Add in a Support to another.
 */
Support& operator+=(Support& one, const Support& other);

/**
 * Scale a Support by a factor.
 */
template<typename Scalar>
Support operator*(const Support& support, const Scalar& scale) {
    Support prod;
    prod.set_forward(support.forward() * scale);
    prod.set_reverse(support.reverse() * scale);
    prod.set_left(support.left() * scale);
    prod.set_right(support.right() * scale);
    
    // log-scaled quality can just be multiplied
    prod.set_quality(support.quality() * scale);
    
    return prod;
}

/**
 * Scale a Support by a factor, in place.
 */
template<typename Scalar>
Support& operator*=(Support& support, const Scalar& scale) {
    support.set_forward(support.forward() * scale);
    support.set_reverse(support.reverse() * scale);
    support.set_left(support.left() * scale);
    support.set_right(support.right() * scale);
    
    // log-scaled quality can just be multiplied
    support.set_quality(support.quality() * scale);
    
    return support;
}

/**
 * Scale a Support by a factor, the other way
 */
template<typename Scalar>
Support operator*(const Scalar& scale, const Support& support) {
    Support prod;
    prod.set_forward(support.forward() * scale);
    prod.set_reverse(support.reverse() * scale);
    prod.set_left(support.left() * scale);
    prod.set_right(support.right() * scale);
    
    // log-scaled quality can just be multiplied
    prod.set_quality(support.quality() * scale);
    
    return prod;
}

/**
 * Divide a Support by a factor.
 */
template<typename Scalar>
Support operator/(const Support& support, const Scalar& scale) {
    Support scaled;
    
    scaled.set_forward(support.forward() / scale);
    scaled.set_reverse(support.reverse() / scale);
    scaled.set_left(support.left() / scale);
    scaled.set_right(support.right() / scale);
    
    // log-scaled quality can just be divided. Maybe.
    scaled.set_quality(support.quality() / scale);
    
    return scaled;
}

/**
 * Divide a Support by a factor, in place.
 */
template<typename Scalar>
Support& operator/=(Support& support, const Scalar& scale) {
    support.set_forward(support.forward() / scale);
    support.set_reverse(support.reverse() / scale);
    support.set_left(support.left() / scale);
    support.set_right(support.right() / scale);
    
    // log-scaled quality can just be divided
    support.set_quality(support.quality() / scale);
    
    return support;
}

/**
 * Support less-than, based on total coverage.
 */
bool operator< (const Support& a, const Support& b);

/**
 * Support greater-than, based on total coverage.
 */
bool operator> (const Support& a, const Support& b);

/**
 * Allow printing a Support.
 */
ostream& operator<<(ostream& stream, const Support& support);

    
}

#endif

