#pragma once

#include <cassert>

#include "DataStructures/Graph/Attributes/AbstractAttribute.h"

// An attribute associating a number of lanes with each edge of a graph.
class NumLanesAttribute : public AbstractAttribute<double> {
 public:
  // Returns the number of lanes of edge e.
  const Type& numLanes(const int e) const {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }

  // Returns a reference to the number of lanes of edge e.
  Type& numLanes(const int e) {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }

 protected:
  static constexpr const char* NAME = "num_lanes"; // The attribute's unique name.
};
