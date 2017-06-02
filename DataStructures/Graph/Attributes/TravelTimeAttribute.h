#pragma once

#include <cassert>

#include "DataStructures/Graph/Attributes/AbstractAttribute.h"
#include "Tools/Constants.h"

// An attribute associating a travel time with each edge of a graph.
class TravelTimeAttribute : public AbstractAttribute<int> {
 public:
  static constexpr Type DEFAULT_VALUE = INFTY;       // The attribute's default value.
  static constexpr const char* NAME = "travel_time"; // The attribute's unique name.

  // Returns the travel time in tenths of a second on edge e.
  const Type& travelTime(const int e) const {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }

  // Returns a reference to the travel time in tenths of a second on edge e.
  Type& travelTime(const int e) {
    assert(e >= 0); assert(e < values.size());
    return values[e];
  }
};
