#pragma once

#include "Tools/Simd/AlignedVector.h"

// The common base class for all attributes. It associates a value of type T with each vertex/edge.
template <typename T>
class AbstractAttribute {
 public:
  using Type = T; // The attribute's underlying type.

  // Returns the attribute's default value.
  static Type defaultValue() {
    return Type();
  }

 protected:
  // Ensures that the attribute can hold at least num values without requiring reallocation.
  void reserve(const int num) {
    values.reserve(num);
  }

  AlignedVector<Type> values; // The attribute's values for the vertices/edges.
};
