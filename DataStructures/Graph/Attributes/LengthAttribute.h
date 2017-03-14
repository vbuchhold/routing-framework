#pragma once

#include <cassert>

#include "DataStructures/Graph/Attributes/AbstractAttribute.h"
#include "Tools/Constants.h"

// An attribute associating a physical length with each edge of a graph.
class LengthAttribute : public AbstractAttribute<int> {
 public:
  static constexpr Type DEFAULT_VALUE = INFTY;    // The attribute's default value.
  static constexpr const char* NAME   = "length"; // The attribute's unique name.

  // Returns the length in meters of edge e.
  const Type& length(const int e) const {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }

  // Returns a reference to the length in meters of edge e.
  Type& length(const int e) {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }
};
