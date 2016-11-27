#pragma once

#include <cassert>

#include "DataStructures/Graph/Attributes/AbstractAttribute.h"

// An attribute associating a physical length with each edge of a graph.
class LengthAttribute : public AbstractAttribute<int> {
 public:
  static constexpr Type DEFAULT_VALUE = INFTY;    // The attribute's default value.
  static constexpr const char* NAME   = "length"; // The attribute's unique name.

  // Returns the length in meters of edge e.
  Type length(const int e) const {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }

  // Sets the length of edge e to the specified value in meters.
  void setLength(const int e, const Type val) {
    assert(e >= 0); assert(e < values.size());
    assert(val >= 0);
    values[e] = val;
  }
};
