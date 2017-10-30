#pragma once

#include <boost/dynamic_bitset.hpp>

// The common base class for all bit attributes. It associates a single bit with each vertex/edge.
class AbstractBitAttribute {
 public:
  using Type = bool; // The attribute's underlying type.

  // Returns the attribute's default value.
  static Type defaultValue() {
    return Type();
  }

 protected:
  // Ensures that the attribute can hold at least num values without requiring reallocation.
  void reserve(const int /*num*/) {
    // boost::dynamic_bitset does not support the operation reserve.
  }

  boost::dynamic_bitset<> values; // The attribute's values for the vertices/edges.
};
