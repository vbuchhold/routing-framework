#pragma once

#include <cassert>

#include "DataStructures/Graph/Attributes/AbstractAttribute.h"

// An attribute associating a free-flow speed (the speed at zero flow) with each edge of a graph.
class FreeFlowSpeedAttribute : public AbstractAttribute<int> {
 public:
  // Returns the free-flow speed in km/h on edge e.
  const Type& freeFlowSpeed(const int e) const {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }

  // Returns a reference to the free-flow speed in km/h on edge e.
  Type& freeFlowSpeed(const int e) {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }

 protected:
  static constexpr const char* NAME = "free_flow_speed"; // The attribute's unique name.
};
