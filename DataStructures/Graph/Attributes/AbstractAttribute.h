#pragma once

#include "Tools/Simd/AlignVector.h"

// The common base class for all attributes. It associates a value of type T with each vertex/edge.
template <typename T>
class AbstractAttribute {
  // Workaround for the GNU compiler. Should not be necessary (and is not for Clang).
  template <typename VertexAttributes,typename EdgeAttributes, bool dynamic>
  friend class Graph;

 public:
  using Type = T; // The attribute's underlying type.

  static constexpr Type DEFAULT_VALUE = Type(); // The attribute's default value.

 protected:
  // Ensures that the attribute can hold at least num values without requiring reallocation.
  void reserve(const int num) {
    values.reserve(num);
  }

  AlignVector<Type> values; // The attribute's values for the vertices/edges.
};
