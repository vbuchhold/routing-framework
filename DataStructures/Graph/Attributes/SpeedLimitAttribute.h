#pragma once

#include <cassert>

#include "DataStructures/Graph/Attributes/AbstractAttribute.h"

// An attribute associating a speed limit with each edge of a graph.
class SpeedLimitAttribute : public AbstractAttribute<int> {
 public:
  // Returns the speed limit in km/h on edge e.
  const Type& speedLimit(const int e) const {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }

  // Returns a reference to the speed limit in km/h on edge e.
  Type& speedLimit(const int e) {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }

 protected:
  static constexpr const char* NAME = "speed_limit"; // The attribute's unique name.
};
