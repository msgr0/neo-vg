#include "genotypekit.hpp"
#include "cactus.hpp"
#include "traversal_finder.hpp"

//#define debug

namespace vg {

using namespace std;



double total(const Support& support) {
    return support.forward() + support.reverse();
}

Support support_min(const Support& a, const Support& b) {
    Support to_return;
    to_return.set_forward(min(a.forward(), b.forward()));
    to_return.set_reverse(min(a.reverse(), b.reverse()));
    to_return.set_quality(min(a.quality(), b.quality()));
    return to_return;
}

Support support_max(const Support& a, const Support& b) {
    Support to_return;
    to_return.set_forward(max(a.forward(), b.forward()));
    to_return.set_reverse(max(a.reverse(), b.reverse()));
    to_return.set_quality(max(a.quality(), b.quality()));    
    return to_return;
}

Support flip(const Support& to_flip) {
    Support flipped = to_flip;
    flipped.set_forward(to_flip.reverse());
    flipped.set_reverse(to_flip.forward());
    return flipped;
}

Support operator+(const Support& one, const Support& other) {
    Support sum;
    sum.set_forward(one.forward() + other.forward());
    sum.set_reverse(one.reverse() + other.reverse());
    sum.set_left(one.left() + other.left());
    sum.set_right(one.right() + other.right());
    
    // log-scaled quality can just be added
    sum.set_quality(one.quality() + other.quality());
    
    return sum;
}

Support& operator+=(Support& one, const Support& other) {
    one.set_forward(one.forward() + other.forward());
    one.set_reverse(one.reverse() + other.reverse());
    one.set_left(one.left() + other.left());
    one.set_right(one.right() + other.right());
    
    // log-scaled quality can just be added
    one.set_quality(one.quality() + other.quality());
    
    return one;
}

bool operator< (const Support& a, const Support& b) {
    return total(a) < total(b);
}

bool operator> (const Support& a, const Support& b) {
    return total(a) > total(b);
}

ostream& operator<<(ostream& stream, const Support& support) {
    return stream << support.forward() << "," << support.reverse();
}


}
