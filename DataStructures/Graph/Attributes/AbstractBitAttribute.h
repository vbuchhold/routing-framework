#pragma once

#include <boost/dynamic_bitset.hpp>

// The common base class for all bit attributes. It associates with each vertex/edge a single bit.
class AbstractBitAttribute {
  // Workaround for the GNU compiler. Should not be necessary (and is not for Clang).
  template <typename VertexAttributes,typename EdgeAttributes, bool dynamic>
  friend class Graph;

 public:
  using Type = bool; // The attribute's underlying type.

  static constexpr Type DEFAULT_VALUE = false; // The attribute's default value.

 protected:
  // Ensures that the attribute can hold num values without requiring reallocation.
  void reserve(const int /*num*/) {
    // boost::dynamic_bitset does not support the operation reserve.
  }

  boost::dynamic_bitset<> values; // The attribute's values for the vertices/edges.
};
