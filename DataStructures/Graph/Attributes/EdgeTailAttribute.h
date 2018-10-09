#pragma once

#include <cassert>

#include "DataStructures/Graph/Attributes/AbstractAttribute.h"
#include "Tools/Constants.h"

// An attribute storing explicitly the tail vertex with each edge of a graph.
class EdgeTailAttribute : public AbstractAttribute<int> {
 public:
  // Returns the attribute's default value.
  static Type defaultValue() {
    return INVALID_VERTEX;
  }

  // Returns the tail vertex of edge e.
  const Type& edgeTail(const int e) const {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }

  // Returns a reference to the tail vertex of edge e.
  Type& edgeTail(const int e) {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }

 protected:
  static constexpr const char* NAME = "edge_tail"; // The attribute's unique name.
};
